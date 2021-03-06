#!/usr/bin/python
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description=
                                     "Filter contigs generated by velvet"
                                     " by their length and/or coverage."
                                     " I use vevelt's naming convention,"
                                     " so I won't work with other assemblers.")
    parser.add_argument('-i', '--input',
                        action='store', type=str, #nargs='+',
                        help='Input file: contigs generated by velvet.')
    parser.add_argument('-c', '--coverage',
                        action='store', type=float, default=0,
                        help='Keep only contigs with coverage above or equal to this number.'
                             ' 0 means no thresholding (default).')
    parser.add_argument('-l', '--length',
                        action='store', type=int, default=0,
                        help='Keep only contigs longer or equal to this number.'
                             ' 0 means no thresholding (default).')
    parser.add_argument('-o', '--output',
                        action='store', type=str, default='filtered_contigs.fa',
                        help='output file (defaults to filtered_contigs.fa)')
    arguments = parser.parse_args()
    return arguments

def check_thresholds(line, length_th, coverage_th):
    line = line.split("_")
    length = int(line[3])
    coverage = float(line[5])
    #print "%i %i wieksze? %s" % (length, length_th, str(length >= length_th))
    return length >= length_th and coverage >= coverage_th

def main():
    arguments = parse_arguments()
    infile = open(arguments.input)
    output = open(arguments.output, 'w')
    save = False
    for line in infile:
        if line.startswith(">"):
            save = check_thresholds(line, arguments.length, arguments.coverage)
        if save:
            output.write(line)
    output.close()

if __name__ == '__main__':
    main()

"""
I assume that contig's name follow this convention:
>NODE_14_length_44_cov_33.840908
>NODE_19_length_1660_cov_21.074097
>NODE_30_length_366_cov_20.945354
>NODE_34_length_5002_cov_20.541782
"""
