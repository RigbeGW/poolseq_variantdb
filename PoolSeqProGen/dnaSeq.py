#!/usr/bin/env python
'''
Module for manipulating DNA sequences.
'''

from Bio.Data import CodonTable
from Bio.Seq import Seq



__author__ = "Rigbe G. Weldatsadik"
__copyright__ = "Copyright (c) 2020 Rigbe G. Weldatsadik"
__license__ = "Apache 2.0"
__version__ = "0.0.1"


def index_genbank_features(gb_record, feature_type, qualifier):
    '''
    copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
    '''

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                # There should only be one protein_id per feature, but there
                # are usually several db_xref entries,
                for value in feature.qualifiers[qualifier]:
                    if value in answer:
                        print("WARNING - Duplicate key %s for %s features \
                        %i and %i") % (value, feature_type, answer[value],
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


def check_stopcodon_index_forward(initial,gbk_fileSeq,geneticCodeID):
    '''
    For CDS on the forward strand it returns the genomic position of the stop codon
    '''
    stop_codon = CodonTable.unambiguous_dna_by_id[geneticCodeID].stop_codons

    # the same as check_stopcodon_index_backward() but here since we are on the
    # forward direction the end of the sequence will be the length of the sequence
    # and we use this to avoid trying to access what is beyond this position

    for i in range(initial, len(gbk_fileSeq), 3):
        if initial + 3 <= len(gbk_fileSeq):
            if str(gbk_fileSeq[i:i + 3]) in stop_codon:
                return i, i + 3
        else:
            return initial, len(gbk_fileSeq)


def check_stopcodon_index_backward(initial,gbk_fileSeq,geneticCodeID):
    '''
    For CDS on the reverse strand, it returns the genomic position of the stop codon
    '''

    # i - 3 >=0 is needed since if this is not true, gbk_fileSeq[i:i-3:-1] will
    # return empty, for example if i(initial) is 2 then the expression
    # gbk_fileSeq[i:i-3:-1] will be gbk_fileSeq[2:-1:-1] which doesn't execute to
    # anything cos in the backward direction we are trying to access beyond the
    # end of the sequence (which in this case is at 0).
    # We have to check all the i in the range not just the initial


    stop_codon_complement = [str(Seq(s).complement()) for s in CodonTable.unambiguous_dna_by_id[geneticCodeID].stop_codons]

    for i in range(initial, 0, -3):
        if i - 3 >= 0:
            # str() is needed since gbk_fileSeq is a seq object
            if str(gbk_fileSeq[i:i - 3:-1]) in stop_codon_complement:
                return i, i - 3
        else:
            return i, 0
