#!/usr/bin/env python
'''
Module for manipulating aligned reads and the variants.
'''

import pysam

__author__ = "Rigbe G. Weldatsadik"
__copyright__ = "Copyright (c) 2020 Rigbe G. Weldatsadik"
__license__ = "Apache 2.0"
__version__ = "0.0.1"



def retrieve_readsAndVars(pos_list,bamFile,chr,combo_allinfo):
    '''
    Retrieves reads from a sorted bam file that span the genomic positions
    in the 'pos_list' argument using pysam view and retrieves the unique combination 
    of variants from the aligned reads.
    '''

    # bam file has to be sorted

    begin_pos = pos_list[0]
    end_pos = pos_list[-1]
    read_pos_alt_list = []
    var_inreads = []
    sam = pysam.view(bamFile, chr + ":" + str(begin_pos) + "-" + str(end_pos))
    for record in sam.strip().split('\n'):
        start_pos = int(record.split('\t')[3])
        flag = int(record.split('\t')[1])
        read = record.split('\t')[9]
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
