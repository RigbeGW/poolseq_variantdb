'''

script for creating the Pool-seq based variant proteome database

'''

import argparse
import csv
import re
import sys
import subprocess
import os
from pyteomics import fasta

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

csv.field_size_limit(sys.maxsize)


def index_genbank_features(gb_record, feature_type, qualifier):
    '''
    
    https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/

    '''

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                # There should only be one protein_id per feature, but there
                # are usually several db_xref entries,
                for value in feature.qualifiers[qualifier]:
                    if value in answer:
                        print "WARNING - Duplicate key %s for %s features \
                        %i and %i" % (value, feature_type, answer[value],
                                      index)
                    else:
                        answer[value] = index

    return answer


def check_if_in_del(cds_var_pos, del_pos_list):

    '''

    Checks if the variant position has been deleted by previous variants and
    adjusts the variant position based on the deleted positions preceeding it

    '''

    if del_pos_list:
        del_pos_list.sort()
        for del_sublist in del_pos_list:
            if cds_var_pos in del_sublist:
                return "True"
            elif cds_var_pos > del_sublist[0]:
                cds_var_pos -= len(del_sublist)
                return cds_var_pos
            else:
                return "False"
    else:
        return "False"


def check_if_in_ins(cds_var_pos, ins_pos_list):

    '''

    checks if there are insert positions before the variant position and
    returns those position lists

    '''

    if ins_pos_list:
        ins_pos_list.sort()
        for ins_sublist in ins_pos_list:
            if cds_var_pos >= ins_sublist[0]:
                cds_var_pos += len(ins_sublist)
                return cds_var_pos
            else:
                return cds_var_pos
    else:
        return cds_var_pos


def check_stopcodon_index_backward(initial):

    '''

    For CDS on the reverse strand, it returns the genomic position of the stop codon

    '''

    # i - 3 >=0 is needed since if this is not true, gbk_fileSeq[i:i-3:-1] will
    # return empty, for example if i(initial) is 2 then the expression
    # gbk_fileSeq[i:i-3:-1] will be gbk_fileSeq[2:-1:-1] which doesn't execute to
    # anything cos in the backward direction we are trying to access beyond the
    # end of the sequence (which in this case is at 0).
    # We have to check all the i in the range not just the initial

    for i in range(initial, 0, -3):
       if i - 3 >= 0:
          # str() is needed since gbk_fileSeq is a seq object
          if str(gbk_fileSeq[i:i - 3:-1]) in stop_codon_complement:
             return i, i - 3
       else:
          return i, 0


def check_stopcodon_index_forward(initial):

    '''

    For CDS on the forward strand it returns the genomic position of the stop codon

    '''

    # the same as check_stopcodon_index_backward() but here since we are on the
    # forward direction the end of the sequence will be the length of the sequence
    # and we use this to avoid trying to access what is beyond this position

    for i in range(initial, len(gbk_fileSeq), 3):
       if initial + 3 <= gbk_fileSeq_len:
          if str(gbk_fileSeq[i:i + 3]) in stop_codon:
             return i, i + 3
       else:
          return initial, gbk_fileSeq_len


def write_fasta(entries):

    '''

    Uses pyteomics's fasta.write() function to write the final fasta output.

    '''

    file = args.fastaFile
    fasta.write(entries, file)


def retrieve_reads(pos_list):

    '''

    Retrieves reads from a sorted bam file that contain the genomic positions
    in the 'pos_list' argument using samtools view and awk programs

    '''

    # bam file has to be sorted

    start_pos = pos_list[0]
    end_pos = pos_list[-1]

    sam = subprocess.Popen(["samtools", "view", args.bamFile, args.genomeID + ":" +
                        str(start_pos) + "-" + str(end_pos)], stdout=subprocess.PIPE)

    awk = subprocess.Popen(["awk", "-v", "OFS=\t", '{print $4,$2,$10}'],
                           stdin=sam.stdout, stdout=subprocess.PIPE)
    sam.stdout.close()  # in order for sam to receive a SIGPIPE if awk exits before sam
    output = awk.communicate()[0]
    return output


def write_WT_protein(k, seq_entry):

    '''

    Writes the sequences of the 'wild type' proteins to the fasta output and
    returns the tryptic peptides of the 'wild type' protein.

    '''

    entries = []

    '''protein_id_cds_index=index_genbank_features(gbk_file,"CDS",
                                                "protein_id")'''

    protein_id_cds_index = index_genbank_features(gbk_file, "CDS",
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
        if args.poolID:
            header = protein_id + '|0_WT|' + str(seq_entry) + '_' + str(args.poolID)
        else:
            header = protein_id + '|0_WT|' + str(seq_entry)
        entries.append((header, immu_tran))
        write_fasta(entries)

        # split by K or R when they are not followd by P
        peptide_list_org_group = re.split(r'([K,R](?!P))', immu_tran)
        peptide_list_org = [peptide_list_org_group[i] + peptide_list_org_group[i + 1]
                         for i in xrange(0, len(peptide_list_org_group) - 2, 2)]

        # since if K or R are the last aa, '' will be included in the list
        if peptide_list_org_group[len(peptide_list_org_group) - 1] != '':
            peptide_list_org.append(
                peptide_list_org_group[len(peptide_list_org_group) - 1])

        return (seq_entry, peptide_list_org)


def write_VariantPep(k, v):

    '''

    Returns the variant peptides to be written to the fasta output

    '''


    '''protein_id_cds_index=index_genbank_features(gbk_file,"CDS",\
                                                "protein_id")'''

    protein_id_cds_index = index_genbank_features(gbk_file, "CDS",
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
                    in_del = check_if_in_del(pos_relativeto_cds, del_range_list)

                    if in_del == "True":
                        modified = False
                        break
                    elif in_del == "False":
                        pos_relativeto_cds = pos_relativeto_cds
                    else:
                        pos_relativeto_cds = in_del

                    pos_relativeto_cds = check_if_in_ins(pos_relativeto_cds, ins_range_list)

                    # adding len(ref) to the last indices in both the insertion
                    # and deletion cases insures we skip and delete as many as
                    # the length of the ref so we avoid repeating the nt in the
                    # insertion case and we delete the right amount of nts in
                    # the deletion case. the len(alt) in the deletion case ensures
                    # we dont delete the trailing nt that are reported in the
                    # vcf but shouldnt be deleted.

                    if(diff > 0):
                        ins_range_list.append(xrange(pos_relativeto_cds +
                                              len(ref), pos_relativeto_cds +
                                              len(alt_allele)))

                    try:
                        # this will create a copy of cds_seq so we dont need to
                        # changs cds_seq to mutable explicitly to modify it
                        # this way
                        mut_cds_seq = mut_cds_seq[:pos_relativeto_cds] + \
                                    alt_allele + \
                                        mut_cds_seq[pos_relativeto_cds + \
                                            len(ref):]
                    except IndexError:
                        modified = False
                        break
                    else:

                        modified = True

                else:  # deletion
                    in_del = check_if_in_del(pos_relativeto_cds, del_range_list)

                    if in_del == "True":
                        modified = False
                        break
                    elif in_del == "False":

                        pos_relativeto_cds = pos_relativeto_cds
                    else:

                        pos_relativeto_cds = in_del

                    pos_relativeto_cds = check_if_in_ins(
                        pos_relativeto_cds, ins_range_list)

                    del_range_list.append(xrange(pos_relativeto_cds + len(alt_allele),
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
                start_index, last_index = check_stopcodon_index_backward(
                    initial)

                # the +1 since from the backward the indices are 1-based rather
                # than 0-based
                lengthmodified_cds_seq = gbk_fileSeq[last_index +
                    1:loc.start] + mut_cds_seq

                lengthmodified_cds_seq = lengthmodified_cds_seq.reverse_complement()

            else:
                cur_pos = len(mut_cds_seq) % 3
                initial = loc.end - cur_pos
                start_index, last_index = check_stopcodon_index_forward(initial)

                lengthmodified_cds_seq = mut_cds_seq + gbk_fileSeq[loc.end:last_index]

            mut_tran = str(lengthmodified_cds_seq.translate(table=11, to_stop=True))
            peptide_list_group = re.split(r'([K,R](?!P))', mut_tran)
            peptide_list = [peptide_list_group[i] + peptide_list_group[i + 1]
                           for i in xrange(0, len(peptide_list_group) - 2, 2)]
            if peptide_list_group[len(peptide_list_group) - 1] != '':
                peptide_list.append(peptide_list_group[len(peptide_list_group) - 1])

            return (protein_id, peptide_list)

        else:
            return None


def retrieve_vars_fromReads(combo_allinfo, reads_file):

    '''

    Retrieves the variants from raw reads.

    '''

    read_pos_alt_list = []
    var_inreads = []
    raw_reads = csv.reader(reads_file.splitlines(), delimiter="\t")
    for row in raw_reads:
        start_pos = int(row[0])
        flag = int(row[1])
        read = row[2]
        for allinfo in combo_allinfo:
            if allinfo[4] == 'snp' and start_pos != 0 and flag < 255:
                try:
                    a = read[abs(allinfo[0] - start_pos)]
                    if a == allinfo[1]:
                        read_pos_alt_list.append(allinfo)
                except IndexError:
                    continue

            elif allinfo[4] == 'ins' and start_pos != 0 and flag < 255:
                try:
                    a = read[abs(allinfo[0] - start_pos):abs(allinfo[0] - start_pos) +
                        allinfo[3] + 1]  # allinfo[3]=len(diff)
                    if a == allinfo[1]:
                        read_pos_alt_list.append(allinfo)
                except IndexError:
                    continue

            elif allinfo[4] == 'del' and start_pos != 0 and flag < 255:
                try:
                    a = read[abs(allinfo[0] - start_pos):abs(allinfo[0] - start_pos) +
                        allinfo[3] + 1]
                    if a == allinfo[6]:
                        read_pos_alt_list.append(allinfo)
                except IndexError:
                    continue

        var_found_inread = sorted(read_pos_alt_list)
        read_pos_alt_list = []

        if var_found_inread and var_found_inread not in var_inreads:
            # when two reads have the same variants so we dont write them twice
            var_inreads.append(var_found_inread)
            yield var_found_inread


if __name__ == '__main__':

    cwd = os.getcwd()
    fasta_file = os.path.join(cwd, 'fasta_output')

    parser = argparse.ArgumentParser(
            description='Create Pool-seq based variant proteome database')

    parser.add_argument('--genomeID', type=str, required=True,
                        help='the ID for the reference genome that is found in the bam')
    parser.add_argument('--genbankFile', type=str, required=True,
                        help='the path to genbank file for the reference genome')
    parser.add_argument('--bamFile', type=str, required=True,
                        help='the path to the sorted alignment bam file for retrieving reads from')
    parser.add_argument('--SnpEffTextOutputFile', type=str, required=True,
                        help='the path to the text output file from SnpEff')
    parser.add_argument('--geneticCodeID', type=int, default=1,
                        help='the genetic code id of your species' \
                        'https: // www.ncbi.nlm.nih.gov / Taxonomy / Utils / wprintgc.cgi')
    parser.add_argument('--fastaFile', type=str, default=fasta_file,
                        help='the path to the fasta output file')
    parser.add_argument('--poolID', type=str, default='',
                        help='Qualifiers to add to the fasta header')

    args = parser.parse_args()

    gbk_file = SeqIO.read(args.genbankFile, 'genbank')
    gbk_fileSeq = gbk_file.seq
    gbk_fileSeq_len = len(gbk_fileSeq)

    #start_codon = CodonTable.unambiguous_dna_by_id[args.geneticCodeID].start_codons
    stop_codon = CodonTable.unambiguous_dna_by_id[args.geneticCodeID].stop_codons
    stop_codon_complement = [str(Seq(s).complement()) for s in stop_codon]
    #stop_codon_rev_comp = [str(Seq(s).reverse_complement()) for s in stop_codon]
    combo = {}
    var_dict = {}
    var_dict_list = []
    aa_vs_var_pos = {}
    aa_pos = []
    pr_withVars = set()

    SnpEff_text = csv.reader(open(args.SnpEffTextOutputFile, "rb"), delimiter='\t')

    next(SnpEff_text, None)
    for row in SnpEff_text:
        var_pos = row[0]
        ref = row[2]
        genotype = row[4]
        alt = row[4]
        protein_id = row[8]  # Actually this is the locus_tag
        if protein_id not in pr_withVars:
            pr_withVars.add(protein_id)
        effect = row[9]

        var_dict[var_pos] = [effect, var_pos, ref, alt]
        var_dict_list.append(var_dict)
        var_dict = {}
        if protein_id not in combo:
            combo[protein_id] = var_dict_list
            var_dict_list = []
        else:
            combo[protein_id].extend(var_dict_list)
            var_dict_list = []

    for key, value in combo.items():
        seq_entry = 0
        peptide_duplicate_tracker = []
        var_pos_index = []
        var_pos_encoding = []
        for_read_comparison_allinfo = []
        for_read_retrieval_pos = []

        try:
            (seq_entry, org_peptide_list) = write_WT_protein(key, seq_entry)
        except TypeError:
            pass

        # each vi is a dict with the var_pos as key and
        # [aa_dict,effect,var_pos,ref,alt,gt] as value
	for vi in value:
	    k = vi.keys()[0]
            v = vi.values()[0]
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
		for_read_comparison_allinfo.append((int(k), v[3], v[2], 2, 'del', v[0],nt_after_del))  

	    else:  # snp
                # (var_pos,alt,ref,len(diff),snp,effect)
		for_read_comparison_allinfo.append((int(k), v[3], v[2], len(v[3]) - len(v[2]),'snp', v[0]))  

	file = retrieve_reads(sorted(for_read_retrieval_pos))
	vars_in_aread = retrieve_vars_fromReads(for_read_comparison_allinfo, file)

        # read_varlist is the list containing variants from a single read
	for read_varlist in vars_in_aread:
            try:
	        (pr_id, peptide_list) = write_VariantPep(key, read_varlist)

	    except TypeError:
		pass

	    else:
		var_pos_index = [idx for idx, val in enumerate(for_read_comparison_allinfo)
                                for s in read_varlist if s == for_read_comparison_allinfo[idx]]

		var_pos_encoding = [1 if i in var_pos_index else 0 for i in
                                   xrange(0, len(for_read_retrieval_pos))]
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
			    elif idx == 0 and len(peptide_list) < 2: # only 1 peptide in peptide_list
				peptide = peptide
                            # can write only the left flanking peptide for the last peptide
			    elif idx == len(peptide_list)-1 :  
				peptide = peptide_list[idx-1]+peptide
			    else:
				peptide = peptide_list[idx-1]+peptide+peptide_list[idx+1]

						
			    seq_entry += 1

                            if args.poolID:
                                header = pr_id + '|' + str(var_pos_encoding_decimal) + '|' + str(seq_entry) + '_' + str(args.poolID)
                            else:
                                header = pr_id + '|' + str(var_pos_encoding_decimal) + '|' + str(seq_entry)

			    entries.append((header,peptide))
			    write_fasta(entries)



    # write protein sequences of proteins that had no variants
                    
    protein_id_cds_index = index_genbank_features(gbk_file,"CDS",\
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
            pr_id=cds_feature.qualifiers['protein_id'][0]
        except KeyError:
            pass
        else:
            if args.poolID:
                fasta_header = pr_id + '|WT_novar' + str(seq_entry) + '_' + str(args.poolID)
            else:
                fasta_header = pr_id + '|WT_novar' + str(seq_entry)
                
            entries.append((fasta_header,immu_tran))
            write_fasta(entries)


   

