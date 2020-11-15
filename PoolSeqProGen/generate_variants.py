#!/usr/bin/env python
'''
Module for creating the Pool-seq based variant proteome database.
It imports the varReads,dnaSeq and pepFasta modules from the PoolSeqProGen package.
'''

import argparse
import csv
import sys
import os

from Bio import SeqIO

from pyteomics import fasta

import pysam

from PoolSeqProGen import varReads
from PoolSeqProGen import dnaSeq
from PoolSeqProGen import pepFasta

csv.field_size_limit(sys.maxsize)


__author__ = "Rigbe G. Weldatsadik"
__copyright__ = "Copyright (c) 2020 Rigbe G. Weldatsadik"
__license__ = "Apache 2.0"
__version__ = "0.0.1"



def main():
    cwd = os.getcwd()
    fasta_file = os.path.join(cwd, 'fasta_output')

    parser = argparse.ArgumentParser(
        description='Create Pool-seq based variant proteome database')

    parser.add_argument('--genbankFile', type=str, required=True,
                        help='the path to genbank file for the reference genome')
    parser.add_argument('--bamFile', type=str, required=True,
                        help='the path to the sorted alignment bam file for retrieving reads from')
    parser.add_argument('--SnpEffTextOutputFile', type=str, required=True,
                        help='the path to the text output file from SnpEff.The chromosome name should be the same as that in the bam file')
    parser.add_argument('--geneticCodeID', type=int, default=1,
                        help='the genetic code id of your species'
                        'https: // www.ncbi.nlm.nih.gov / Taxonomy / Utils / wprintgc.cgi')
    parser.add_argument('--fastaFile', type=str, default=fasta_file,
                        help='the path to the fasta output file')
    parser.add_argument('--poolID', type=str, default='',
                        help='Qualifiers to add to the fasta header')

    args = parser.parse_args()

    gbk_file = SeqIO.read(args.genbankFile, 'genbank')
    gbk_fileSeq = gbk_file.seq

    combo = {}
    var_dict = {}
    var_dict_list = []
    aa_vs_var_pos = {}
    aa_pos = []
    pr_withVars = set()

    with open(args.SnpEffTextOutputFile, "r") as snpf:
        SnpEff_text = csv.DictReader(snpf, delimiter='\t')

        for row in SnpEff_text:
            var_pos = row['POS']
            ref = row['REF']
            genotype = row['EFF[*].GT']
            alt = row['EFF[*].GT']
            protein_id = row['EFF[*].TRID']  
            if protein_id not in pr_withVars:
                pr_withVars.add(protein_id)
            effect = row['EFF[*].EFFECT']

            var_dict[var_pos] = [effect, var_pos, ref, alt]
            var_dict_list.append(var_dict)
            var_dict = {}
            if protein_id not in combo:
                combo[protein_id] = var_dict_list
                var_dict_list = []
            else:
                combo[protein_id].extend(var_dict_list)
                var_dict_list = []

        chrom = row['CHROM']

    for key, value in combo.items():
        seq_entry = 0
        peptide_duplicate_tracker = []
        var_pos_index = []
        var_pos_encoding = []
        for_read_comparison_allinfo = []
        for_read_retrieval_pos = []

        try:
            (seq_entry, org_peptide_list) = pepFasta.WriteWtPro_ReturnWtPep(key,seq_entry,gbk_file,args.fastaFile,args.poolID)
        except TypeError:
            pass

        # each vi is a dict with the var_pos as key and
        # [aa_dict,effect,var_pos,ref,alt,gt] as value
        for vi in value:
            k = list(vi.keys())[0]
            v = list(vi.values())[0]
            # holds all the var_pos associated with the protein
            for_read_retrieval_pos.append(int(k))

            if len(v[3]) > len(v[2]):  # ins
                # (var_pos,alt,ref,len(diff),ins,effect)
                for_read_comparison_allinfo.append((int(k), v[3], v[2],
                                                    len(v[3]) - len(v[2]), 'ins', v[0]))

            elif len(v[3]) < len(v[2]):  # del
                # chose 2(i.e. will investigate only 2 nts after the deleted part)
                # to make sure not to miss reads that contain the del and then
                # nearby snps or other indels

                nt_after_del = v[3] + str(gbk_fileSeq[(int(k) - 1) + len(v[2]):int(k) - 1 +
                                                      len(v[2]) + 2])
                # (var_pos,alt,ref,2,del,effect,nt_after_del)
                for_read_comparison_allinfo.append(
                    (int(k), v[3], v[2], 2, 'del', v[0], nt_after_del))

            else:  # snp
                # (var_pos,alt,ref,len(diff),snp,effect)
                for_read_comparison_allinfo.append(
                    (int(k), v[3], v[2], len(v[3]) - len(v[2]), 'snp', v[0]))
        vars_in_aread = varReads.retrieve_readsAndVars(sorted(for_read_retrieval_pos),args.bamFile,chrom,
            for_read_comparison_allinfo)

        # read_varlist is the list containing variants from a single read
        for read_varlist in vars_in_aread:
            try:
                (pr_id, peptide_list) = pepFasta.ReturnVarPep(key,read_varlist,gbk_file,gbk_fileSeq,args.geneticCodeID)

            except TypeError:
                pass

            else:
                var_pos_index = [idx for idx, val in enumerate(for_read_comparison_allinfo)
                                 for s in read_varlist if s == for_read_comparison_allinfo[idx]]

                var_pos_encoding = [1 if i in var_pos_index else 0 for i in
                                    range(0, len(for_read_retrieval_pos))]
                # bitwise encoding
                var_pos_encoding_bitwise = ''.join(map(str, var_pos_encoding))
                # the bitwise encoding changed to decimal representation for brevity
                var_pos_encoding_decimal = int(var_pos_encoding_bitwise, 2)

                if peptide_list is not None:
                    for idx, peptide in enumerate(peptide_list):
                        entries = []
                        if peptide not in peptide_duplicate_tracker and peptide not in org_peptide_list and len(peptide) >= 4:
                            peptide_duplicate_tracker.append(peptide)
                            # can write only the right flanking peptide for the first peptide
                            if idx == 0 and len(peptide_list) >= 2:
                                peptide += peptide_list[idx+1]
                            # only 1 peptide in peptide_list
                            elif idx == 0 and len(peptide_list) < 2:
                                peptide = peptide
                            # can write only the left flanking peptide for the last peptide
                            elif idx == len(peptide_list)-1:
                                peptide = peptide_list[idx-1]+peptide
                            else:
                                peptide = peptide_list[idx-1] + \
                                    peptide+peptide_list[idx+1]

                            seq_entry += 1

                            if args.poolID:
                                header = pr_id + '|' + \
                                    str(var_pos_encoding_decimal) + '|' + \
                                    str(seq_entry) + '_' + str(args.poolID)
                            else:
                                header = pr_id + '|' + \
                                    str(var_pos_encoding_decimal) + \
                                    '|' + str(seq_entry)

                            entries.append((header, peptide))
                            fasta.write(entries, args.fastaFile)

                            

    # write protein sequences of proteins that had no variants

    protein_id_cds_index = dnaSeq.index_genbank_features(gbk_file, "CDS",
                                                  "locus_tag")
    original_pr = set(protein_id_cds_index.keys())

    novar_pr = original_pr.difference(pr_withVars)

    for pr in novar_pr:
        entries = []
        seq_entry += 1
        index = protein_id_cds_index[pr]
        cds_feature = gbk_file.features[index]
        try:
            immu_tran = cds_feature.qualifiers['translation'][0]
            pr_id = cds_feature.qualifiers['protein_id'][0]
        except KeyError:
            pass
        else:
            if args.poolID:
                fasta_header = pr_id + '|WT_novar' + \
                    str(seq_entry) + '_' + str(args.poolID)
            else:
                fasta_header = pr_id + '|WT_novar' + str(seq_entry)

            entries.append((fasta_header, immu_tran))
            fasta.write(entries, args.fastaFile)
            



main()