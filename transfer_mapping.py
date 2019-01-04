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
        #cigar = contig.cigartuples
        start = contig.query_alignment_start
        end = start + contig.query_alignment_length
        chromosome = contig.reference_name
        start_ref, end_ref = contig.reference_start, contig.reference_end
        is_reverse = contig.is_reverse
        dictionary[contig][(start, end)] = (chromosome, start_ref, end_ref,
                                            is_reverse)
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
        read.reference_start = contig_end + start

def find_contig_mapping(mapping, start=0):
    for on_contig, on_reference in mapping.items():
        if on_contig[0] >= start and on_contig[1] < start:
            return on_reference
    return None

#def daj_nazwe_readu(read):
#    return read.qname + str(int(read.is_read1))

def merge_files(reads2reference, reads2contigs, contigs2reference,
                outfile):
    new_bam = pysam.AlignmentFile(outfile, "wb")
    for nr, bam in enumerate(reads2contigs):
        contig_mapping = read_in_contig_mapping(contigs2reference[nr])
        for read in bam:
            transfer_read(read, contig_mapping)
            new_bam.write(read)
    for read in reads2reference:
        new_bam.write(read)

def main():
    arguments = parse_arguments()
    reads2reference = pysam.AlignmentFile(arguments.reads2reference)
    reads2contigs = [pysam.AlignmentFile(filename) for filename in arguments.reads2contigs]
    contigs2reference = [pysam.AlignmentFile(filename) for filename in arguments.showcoords]
    merge_files(reads2reference, reads2contigs, contigs2reference,
                arguments.output)

if __name__ == '__main__':
    main()
