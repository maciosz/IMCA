#!/usr/bin/python
import argparse
import sys
import os
import collections
import pysam


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='reads2reference', action = 'store', type=str,
                       help='.bam / .sam mapped to reference')
    parser.add_argument('-c', dest='reads2contigs', action = 'store', type=str, nargs='+',
                        help = '.bam / .sam mapped to contigs')
    parser.add_argument('-m', dest='contigs2reference', action='store', type=str, nargs='+',
                       help='file with contigs mapped to reference')
    parser.add_argument('-d', dest='output_dir', action='store', type=str, default='merged',
    			help = 'destination to save output files (defaults to merged)')
    parser.add_argument('--suf', dest='output_sufix', action='store', type=str, default='',
                        help='suffix to output files (defaults to None)')
    arguments = parser.parse_args()
    check_arguments(arguments)
    return arguments

def check_arguments(arguments):
    pass

"""
Zakladam ze bede miala z gtfa / sama contigs2reference taki slownik:
{contig: opis contigu}
gdzie opis contigu to slownik:
{obszar: obszar_na_ref}
Np.:
{'contigA':
  {(1, 10): ('chr1', 40, 50),
   (10, 50): ('chr1', 100, 140)},
 'contigB': None,
 'contigC':
  {(1, 1000): ('chr3', 4000, 4999),
   (1000, 2000): ('chr5', 1, 999)}
}
"""

def transfer_read(read):
    contig, start, end = read.reference_name, read.reference_start, read.reference_end
    contig_coords = find_contig_mapping(contig, start)
    chromosome, contig_start, contig_end = contig_coords
    read.reference_start = contig_start + start
    read.reference_id = chromosome
    # powinnam uwzgledniac ze contig mogl sie zmapowac odwrotnie

def find_contig_mapping(contig, start=0):
    mapping = contigs2reference[contig]
    for on_contig, on_reference in mapping.items():
        if on_contig[0] >= start and on_contig[1] < start:
            return on_reference
    return None

def read_in_contig_mapping(infile):
    return

#def daj_nazwe_readu(read):
#    return read.qname + str(int(read.is_read1))

def merge_files():
    new_bam = pysam.AlignmentFile(outfile, "wb")
    for bam in infiles:
        for read in bam:
            transfer_read(read)
            new_bam.write(read)


def zmerguj_zmapowane_pliki(reads2reference, reads2contigs,
                        zmapowanie_contigow, showsnps,
                        output_dir, output_sufix):
    merged = pysam.AlignmentFile(output_dir + "/merged" + output_sufix + ".bam", "wb", template=reads2reference)
    zapisane = set()
    for read in reads2contigs:
        if read.is_unmapped:
            continue
        zmapowanie = przenies_mapowanie(read, zmapowanie_contigow, showsnps)
        if zmapowanie:
            merged.write(read)
            zapisane.add(daj_nazwe_readu(read))
    for read in reads2reference:
        if read not in zapisane:
            merged.write(read)
    merged.close()

def main():
    arguments = parse()
    reads2reference = pysam.AlignmentFile(arguments.reads2reference)
    reads2contigs = [pysam.AlignmentFile(filename) for filename in arguments.reads2contigs][0]
    zmapowanie_contigow = [wczytaj_showcoords(filename) for filename in arguments.showcoords][0]
    showsnps = [wczytaj_showsnps(filename) for filename in arguments.showsnps][0]
    if not os.path.isdir(arguments.output_dir):
        os.mkdir(argumenty.output_dir)
    zmerguj_zmapowane_pliki(reads2reference, reads2contigs,
                        zmapowanie_contigow, showsnps, 
                        argumenty.output_dir, argumenty.output_sufix)
 

def write_merged_file(mapped2original, reads2contigs, showcoords, output_dir, output_suffix):
    merged_bam = pysam.AlignmentFile(output_dir + '/merged' + output_suffix + '.bam', 'wb', template=mapped2original)
    saved_reads = []
    for numer_pliku in xrange(len(reads2contigs)):
        print "Transfering mapping of reads for file %d" % numer_pliku
        przenies_mapowanie(reads2contigs[numer_pliku], showcoords[numer_pliku], merged_bam, saved_reads)
    print "zapisuje ready zmapowane na original"
    zapisz_ready(mapped2original, merged_bam, saved_reads)
    merged_bam.close()

def przenies_mapowanie(bam, showcoords, output, zapisane_ready):
    ile_niezmienionych = 0
    ile_zmienionych = 0
    for read in bam:
        if daj_nazwe_readu(read) in zapisane_ready:
            continue
        czy_zmienione = zmien_koordynaty(read, showcoords)
        if not czy_zmienione:
            ile_niezmienionych += 1
            if ile_niezmienionych % 1000 == 0:
                print "niezmienionych:", ile_niezmienionych
        if czy_zmienione:
            ile_zmienionych += 1
            if ile_zmienionych % 10000 == 0:
                print "zmienionych:", ile_zmienionych
            output.write(read)
            zapisane_ready.append(daj_nazwe_readu(read))
    print "nie zmienionych:", ile_niezmienionych
    
def zmien_koordynaty(read, showcoords):
    start_read2contig = read.reference_start
    end_read2contig = read.reference_end
    is_reversed = False
    znalazlem = False
    for linia in showcoords:
        start_ref, end_ref, start_contig, end_contig = linia
        if start_contig > end_contig:
            is_reversed = True
            start_contig, end_contig = end_contig, start_contig
        if start_read2contig >= start_contig and start_read2contig <= end_contig:
            znalazlem = True
            break
    if not is_reversed:
        start_read2original = start_read2contig + start_ref
    elif is_reversed:
        start_read2original = start_ref + (start_contig - end_read2contig)
    #if not znalazlem:

if __name__ == '__main__':
    main()



