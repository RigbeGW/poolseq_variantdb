#!/usr/bin/env python
'''
Module for in silico digestion of WT and Variant proteins and writing to fasta files.
'''

import re

from pyteomics import fasta

from PoolSeqProGen import dnaSeq


__author__ = "Rigbe G. Weldatsadik"
__copyright__ = "Copyright (c) 2020 Rigbe G. Weldatsadik"
__license__ = "Apache 2.0"
__version__ = "0.0.1"



def WriteWtPro_ReturnWtPep(k, seq_entry,gbk_file,fastaOutput,poolID):
    '''
    Writes the sequences of the 'wild type' proteins to the fasta output and
    returns the tryptic peptides of the 'wild type' protein.
    '''

    entries = []

    '''protein_id_cds_index=dnaSeq.index_genbank_features(gbk_file,"CDS",
                                                "protein_id")'''

    protein_id_cds_index = dnaSeq.index_genbank_features(gbk_file, "CDS",
                                                  "locus_tag")

    index = protein_id_cds_index[k]
    cds_feature = gbk_file.features[index]
    try:  # there are cds with no translation and protein_id
        immu_tran = cds_feature.qualifiers['translation'][0]
        protein_id = cds_feature.qualifiers['protein_id'][0]
    except KeyError:
        pass
    else:
        seq_entry += 1
        if poolID:
            header = protein_id + '|0_WT|' + \
                str(seq_entry) + '_' + str(poolID)
        else:
            header = protein_id + '|0_WT|' + str(seq_entry)

        entries.append((header, immu_tran))
        fasta.write(entries, fastaOutput)

        
        # split by K or R when they are not followd by P
        peptide_list_org = re.compile(".(?:(?<![KR](?!P)).)*").findall(str(immu_tran))


        return (seq_entry, peptide_list_org)


def ReturnVarPep(k,v,gbk_file,gbk_fileSeq,geneticCodeID):
    '''
    Returns the variant peptides to be written to the fasta output
    '''

    '''protein_id_cds_index=dnaSeq.index_genbank_features(gbk_file,"CDS",\
                                                "protein_id")'''

    protein_id_cds_index = dnaSeq.index_genbank_features(gbk_file, "CDS",
                                                  "locus_tag")

    index = protein_id_cds_index[k]
    cds_feature = gbk_file.features[index]
    try:
        immu_tran = cds_feature.qualifiers['translation'][0]
        protein_id = cds_feature.qualifiers['protein_id'][0]
    except KeyError:
        pass
    else:
        loc = cds_feature.location
        strand = cds_feature.location.strand  # +1 or -1
        feature_seq = cds_feature.extract(gbk_file.seq)
        pr_name = cds_feature.qualifiers['product'][0]

        # this wont return the reverse complement for protein on the complement strand
        # unlike extract method
        mut_cds_seq = gbk_file[loc.start:loc.end].seq
        mut_cds_len = len(mut_cds_seq)
        del_range_list = []
        ins_range_list = []
        modified = False

        # subset is in the form of (var_pos,alt,ref,effect,snp) for snp
        # (var_pos,alt,ref,effect,ins) for insertion
        # (var_pos,nt_after_del,alt,ref,del,effect) for deletion
        for subset in v:
            var_pos = subset[0]
            alt_allele = subset[1]
            ref = subset[2]
            effect = subset[5]

            # since var_pos is 1-based position from SnpEff and relative to
            # the whole genome length
            zero_based_pos = int(var_pos) - 1
            pos_relativeto_cds = abs(loc.start - zero_based_pos)

            diff = len(alt_allele) - len(ref)

            if effect == "START_LOST":
                modified = False
                break

            if diff >= 0:  # ins and snp
                # checking deletion range list first since if that particular
                # position has been deleted, it doesn't make sense to
                # insert there
                in_del = dnaSeq.check_if_in_del(pos_relativeto_cds, del_range_list)

                if in_del == "True":
                    modified = False
                    break
                elif in_del == "False":
                    pos_relativeto_cds = pos_relativeto_cds
                else:
                    pos_relativeto_cds = in_del

                pos_relativeto_cds = dnaSeq.check_if_in_ins(
                    pos_relativeto_cds, ins_range_list)

                # adding len(ref) to the last indices in both the insertion
                # and deletion cases insures we skip and delete as many as
                # the length of the ref so we avoid repeating the nt in the
                # insertion case and we delete the right amount of nts in
                # the deletion case. the len(alt) in the deletion case ensures
                # we dont delete the trailing nt that are reported in the
                # vcf but shouldnt be deleted.

                if(diff > 0):
                    ins_range_list.append(range(pos_relativeto_cds +
                                                 len(ref), pos_relativeto_cds +
                                                 len(alt_allele)))

                try:
                    # this will create a copy of cds_seq so we dont need to
                    # changs cds_seq to mutable explicitly to modify it
                    # this way
                    mut_cds_seq = mut_cds_seq[:pos_relativeto_cds] + \
                        alt_allele + \
                        mut_cds_seq[pos_relativeto_cds +
                                    len(ref):]
                except IndexError:
                    modified = False
                    break
                else:

                    modified = True

            else:  # deletion
                in_del = dnaSeq.check_if_in_del(pos_relativeto_cds, del_range_list)

                if in_del == "True":
                    modified = False
                    break
                elif in_del == "False":

                    pos_relativeto_cds = pos_relativeto_cds
                else:

                    pos_relativeto_cds = in_del

                pos_relativeto_cds = dnaSeq.check_if_in_ins(
                    pos_relativeto_cds, ins_range_list)

                del_range_list.append(range(pos_relativeto_cds + len(alt_allele),
                                             pos_relativeto_cds + len(ref)))

                try:
                    mut_cds_seq = mut_cds_seq[:pos_relativeto_cds + len(alt_allele)] + \
                        mut_cds_seq[pos_relativeto_cds + len(ref):]
                except IndexError:
                    modified = False
                    break
                else:

                    modified = True

        if modified:
            if strand == -1:
                cur_pos = len(mut_cds_seq) % 3
                if cur_pos == 0:
                    initial = loc.start - 1
                elif cur_pos == 1:
                    initial = loc.start  # since there is a nt at loc.start of the cds
                else:
                    initial = loc.start + 1
                start_index, last_index = dnaSeq.check_stopcodon_index_backward(
                    initial,gbk_fileSeq,geneticCodeID)

                # the +1 since from the backward the indices are 1-based rather
                # than 0-based
                lengthmodified_cds_seq = gbk_fileSeq[last_index +
                                                     1:loc.start] + mut_cds_seq

                lengthmodified_cds_seq = lengthmodified_cds_seq.reverse_complement()

            else:
                cur_pos = len(mut_cds_seq) % 3
                initial = loc.end - cur_pos
                start_index, last_index = dnaSeq.check_stopcodon_index_forward(
                    initial,gbk_fileSeq,geneticCodeID)

                lengthmodified_cds_seq = mut_cds_seq + gbk_fileSeq[loc.end:last_index]

            mut_tran = str(lengthmodified_cds_seq.translate(
                table=11, to_stop=True))
            peptide_list = re.compile(".(?:(?<![KR](?!P)).)*").findall(mut_tran)

            return (protein_id, peptide_list)

        else:
            return None

        

