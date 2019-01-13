#!/usr/bin/python
#import sys
#import os
import argparse
import collections
import pysam


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads2reference',
                        action='store', type=str,
                        help='.bam / .sam mapped to reference')
    parser.add_argument('-c', '--reads2contigs',
                        action='store', type=str, nargs='+',
                        help='.bam / .sam mapped to contigs')
    parser.add_argument('-m', '--contigs2reference',
                        action='store', type=str, nargs='+',
                        help='bam / sam with contigs mapped to reference')
    parser.add_argument('-o', '--output',
                        action='store', type=str, default='merged.bam',
                        help='output file (defaults to merged.bam)')
    parser.add_argument('-t', '--transfer-from-reference',
                        action='store_true', #type=bool, default=False, const=True,
                        help='transfer reads mapped to the reference?'
                        ' By default I won\'t.')
    parser.add_argument('-k', '--keep-unmapped-contigs',
                        action='store_true', # type=bool, default=False, const=True,
                        help='Add unmapped contigs to the references in the output?'
                        ' By default I won\'t.')
    parser.add_argument('-a', '--ambigous-mappings',
                        action='store', type=str, default='best',
                        help='What to do with ambigous mappings?\n'
                        ' - best: keep the best mapping (default).'
                        '         if there is none, choose randomly.\n\n'
                        ' - best-rm: keep the best mapping.'
                        '         if there is none, remove this contig/read.\n'
                        '- rm: remove the ambigously mapped reads/contigs.\n'
                        '- keep: keep all the best mappings.\n')
    arguments = parser.parse_args()
    return arguments

#{contig: opis contigu}
#gdzie opis contigu to slownik:
#{obszar: obszar_na_ref}
#Np.:
#{'contigA':
#  {(1, 10): ('chr1', 40, 50),
#   (10, 50): ('chr1', 100, 140)},
# 'contigB': None,
# 'contigC':
#  {(1, 1000): ('chr3', 4000, 4999),
#   (1000, 2000): ('chr5', 1, 999)}
#}

def read_in_contig_mapping(filename):
    bam = pysam.AlignmentFile(filename)
    dictionary = collections.defaultdict(dict)
    for contig in bam:
        if contig.is_unmapped:
            dictionary[contig] = None
            continue
        cigar = contig.cigartuples
        start = contig.query_alignment_start
        end = start + contig.query_alignment_length
        chromosome = contig.reference_name
        start_ref, end_ref = contig.reference_start, contig.reference_end
        is_reverse = contig.is_reverse
        previous_start = start
        previous_end = start
        for what, how_many in cigar:
            if what == 0 or what == 7 or what == 8: # match / mismatch
                block_start = previous_end
                block_end = block_start + how_many
                previous_start, previous_end = block_start, block_end
                dictionary[contig][(block_start, block_end)] = (chromosome, start_ref, end_ref,
                                                                is_reverse)
            elif what == 1: # insertion
                block_start = previous_start
                block_end = previous_end
                dictionary[contig][(block_start, block_end)] = (chromosome, start_ref, end_ref,
                                                                is_reverse)
            else:   # deletion
                previous_start = previous_end
                previous_end += how_many
    return dictionary

def transfer_read(read, contig_mapping):
    contig, start = read.reference_name, read.reference_start
    contig_mapping = contig_mapping[contig]
    contig_coords = find_contig_mapping(contig_mapping, start)
    if contig_coords is None:
        return # ?
    chromosome, contig_start, contig_end, is_reverse = contig_coords
    read.reference_id = chromosome
    if not is_reverse:
        read.reference_start = contig_start + start
    else:
        read.reference_start = contig_end - start
        read.is_reverse = not read.is_reverse

def find_contig_mapping(mapping, start=0):
    for on_contig, on_reference in mapping.items():
        if on_contig[0] >= start and on_contig[1] < start:
            return on_reference
    return None

def merge_files(reads2reference, reads2contigs, contigs2reference,
                arguments):
    outfile = arguments.output
    new_bam = pysam.AlignmentFile(outfile, "wb", template = reads2reference)
    for nr, bam in enumerate(reads2contigs):
        contig_mapping = read_in_contig_mapping(contigs2reference[nr])
        for read in bam:
            transfer_read(read, contig_mapping)
            new_bam.write(read)
    if arguments.transfer_from_reference:
        for nr, bam in enumerate(reads2contigs):
            for read in reads2reference:
                transfer_read(read, contig_mapping)
                new_bam.write(read)
    else:
        for read in reads2reference:
            new_bam.write(read)

def main():
    arguments = parse_arguments()
    reads2reference = pysam.AlignmentFile(arguments.reads2reference)
    reads2contigs = [pysam.AlignmentFile(filename) for filename in arguments.reads2contigs]
    contigs2reference = [pysam.AlignmentFile(filename) for filename in arguments.showcoords]
    merge_files(reads2reference, reads2contigs, contigs2reference,
                arguments)

if __name__ == '__main__':
    main()
