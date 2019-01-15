# IMCA

## Improving Mapping with Contigs Assembly

When mapping reads from any NGS experiment
 (like ChIP-seq, DNase-seq etc.)
 there is always some loss due to the differences
 between the reference genome and the actual genome being sequenced.
 IMCA attempts to reduce this loss.
 In short, first contigs are assembled from sequenced reads -
 they should correspond to the actual sequence of the genome.
 Next, the contigs are aligned to the reference genome,
 to find out where those sequences come from.
 Then reads are mapped to them.
 If there are any reads that didn't map to the reference
 due to differences in sequence,
 it might be that they will map to some contig
 that will map to the reference;
 then we can transfer coordinates of this read's mapping
 to the reference.
 (picture here would be nice)

If you encounter any problems running IMCA
 or have any comments or suggestions,
 please let me know at a dot macioszek at mimuw.edu.pl.

### Instalation and requirements

No installation needed.
At least not for these scripts.
You might need to add executive rights using `chmod +x` command to the python scripts.
Also you might need to install some software needed to complete additional steps.
Specifically,
you need some mapper for your reads, mapper for contigs, and assembler to create contigs.
The script was tested using minimap2 and velvet,
but should work with output of any other software,
as long as formats are correct:
the mappings should be in .bam / .sam format
(preferalby with CIGAR strings assigned to every mapping),
and contigs in fasta format.

`transfer_mapping.py` script requires python modules argparse and pysam.
 It was tested with versions argparse1.2.1 and pysam0.9.1.4,
 python2.7.6.
`filter_velvet_contigs.py` needs argparse.

### Workflow

1. Map reads to the reference.
2. Construct contigs from the reads.
3. Potentially, filter contigs.
4. Map (align) contigs to the reference.
5. Map reads to the mapped contigs.
6. Transfer the mapping of the reads from contigs to reference.

Steps from 1 to 5 can be done with any software of your choice.
Step 6 can be done using `transfer_mapping.py` script.
Below I present example workflow.

1. Map reads to the reference.

Using [minimap2](https://github.com/lh3/minimap2):

`minimap2 -a -x sr reference.fa reads.fq > reads2reference.sam 2> mapping.err`

You can also use any other mapper, like bowtie2, BWA, Tophat etc.

2. Construct contigs from the reads.

Using [velvet](https://www.ebi.ac.uk/~zerbino/velvet/):

`velveth velvet_output 27 -fastq reads.fq`
`velvetg velvet_output`

This will generate directory `velvet_output` with several files,
including large Sequences, Roadmap, Graph and LastGraph, that you can remove
once `velvetg` ends running.
We are only going to need contigs.fa and potentially stats.txt for filtering.
`27` refers to the size of the k-mers constructed by velvet
and of course can be changed according to your needs.
Please refer to [velvet's documantation](https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf)
for more details.

You can use any other assembler.
The important thing is to get contigs in fasta format in this step.

3. Potentially, filter contigs.

You might choose to filter out contigs that are too short,
or have too small coverage.
I encourage to look at the obtained contigs and decide what seems best
to obtain better results.

If you used velvet as assembler, you may use `filter_velvet_contigs.py` script.
It doesn't work with other assemblers, because it uses
velvet's naming convention.
It can filter by contig's length and coverage.
For example, run:

`./filter_velvet_contigs.py -i velvet/contigs.fa -o velvet/filtered_contigs.fa -l 500 -c 5`

to create file `filtered_contigs.fa`, that will contain
these contigs from `velvet/contigs.fa` that were at least 500 bp long
and had coverage higher or equal to 5.

4. Map (align) contigs to the reference.

Using minimap2:

`minimap2 -a -x asm10 reference.fa velvet_output/contigs.fa > contigs2reference.sam 2> mapping_contigs.err`

You can also use any other mapper capable of mapping large contigs
or tool designed specifically for aligning genomes,
like MUMmer, BLAST etc.

Important note: the output of this step must be in sam or bam format.
If you use a software that generates anything else,
you need to convert it to sam / bam format before last step.
(For MUMmer, check out this [link](https://www.biostars.org/p/185384/).)
For better accuracy, CIGAR strings should be available for every alignment.
See documentation for details.

5. Map reads to the mapped contigs.

This task is simmilar to the step 1; we just use mapped and potentially filtered contigs instead of reference.
For example, using minimap2:

`minimap2 -a -x sr velvet_output/contigs.fa reads.fq > reads2contigs.sam 2> mapping2contigs.err`

But you can use Bowtie2, BWA, Tophat etc.

6. Transfer the mapping of the reads from contigs to reference.

`./transferMapping.py -r reads2reference.sam -c reads2contigs.sam -m contigs2reference.sam -o merged_mapping.bam`

This will generate file `merged_mapping.bam`
 with all the reads that were mapped to the reference
 and reads mapped to the contigs which mapped to the reference.
 By default, the reads mapped to the reference stay unchanged,
 but you may choose to change their coordinates when possible by adding `-t` argument.
 See documentation for more details.

### Parameters

To see all the arguments available, run the script with `-h` option:

```
usage: transfer_mapping.py [-h] [-r READS2REFERENCE]
                           [-c READS2CONTIGS [READS2CONTIGS ...]]
                           [-m CONTIGS2REFERENCE [CONTIGS2REFERENCE ...]]
                           [-o OUTPUT] [-t] [-k] [-a AMBIGOUS_MAPPINGS]

optional arguments:
  -h, --help            show this help message and exit
  -r READS2REFERENCE, --reads2reference READS2REFERENCE
                        .bam / .sam mapped to reference
  -c READS2CONTIGS [READS2CONTIGS ...], --reads2contigs READS2CONTIGS [READS2CONTIGS ...]
                        .bam / .sam mapped to contigs
  -m CONTIGS2REFERENCE [CONTIGS2REFERENCE ...], --contigs2reference CONTIGS2REFERENCE [CONTIGS2REFERENCE ...]
                        bam / sam with contigs mapped to reference
  -o OUTPUT, --output OUTPUT
                        output file (defaults to merged.bam)
  -t, --transfer-from-reference
                        transfer reads mapped to the reference? By default I
                        won't.
  -k, --keep-unmapped-contigs
                        Add unmapped contigs to the references in the output?
                        By default I won't.
  -a AMBIGOUS_MAPPINGS, --ambigous-mappings AMBIGOUS_MAPPINGS
                        What to do with ambigous mappings? - best: keep the
                        best mapping (default). if there is none, choose
                        randomly. - best-rm: keep the best mapping. if there
                        is none, remove this contig/read. - rm: remove the
                        ambigously mapped reads/contigs. - keep: keep all the
                        best mappings.
```

#### Mandatory arguments:

`-r` or `--reads2reference` is the file with reads mapped to the reference.
 It should be in bam or sam format and the file should have index generated;
 you can do it using `samtools index` command.

`-c` or `--reads2contigs` is the file with reads mapped to the contigs.
 It should be in bam or sam format and the file should have index generated;
 you can do it using `samtools index` command.
 To ensure better estimation of transfer coordinates,
 every mapped read should have CIGAR string set.
 It is not compulsory, but it is highly encouraged.
 Many mappers will generate them by default
 (Bowtie2, minimap2, BWA...).

`-m` or `--contigs2reference` is the file with contigs mapped to the reference
 (`m` stands for mapping, because `c` and `r` were already taken...).
 It should be in bam or sam format and the file should have index generated;
 you can do it using `samtools index` command.
 To ensure better estimation of transfer coordinates,
 every mapped contig should have CIGAR string set.
 It is not compulsory, but it is highly encouraged.
 If you use some tool designed for aligning contigs to the reference
 it might need some additional parsing to retrieve CIGAR strings
 from it's output.
 You can read more about how CIGAR string is constructed [here](https://www.drive5.com/usearch/manual/cigar.html).

#### Optional arguments:

`-o` is the desired named of the output file. Defaults to `merged.bam`.
 Currently output is always in .bam format.

`-t` or `--transfer-from-reference` is a flag that decides
 what to do with reads that are mapped both to the contigs and to the reference:
 should IMCA leave them alone and keep the coordinates from `reads2reference` file
 (that's the default behaviour)
 or should their coordinates be transferred when possible?
 The latter happens if you run the script with this option.

`-k` or `--keep-unmapped-contigs` is a flag that decides
 what to do with reads that are mapped to contigs which didn't map to the reference:
 should we just ignore and remove them
 (that's the default behaviour)
 or should we add these contigs as references to the output .bam file
 and keep the reads?
 The latter happens if you run the script with this option.
 Keep in mind that the header of your output bam file
 will now be different than the header of `reads2reference`;
 it will contain some additional references.

`-a` or `--ambigous-mappings` decides
 what to do with contigs and reads that are mapped to multiple positions.
 This argument requires a value describing the desired behaviour.
 `best` means: keep the best mapping;
 if there are multiple equally good mappings,
 randomly choose one (that's the default).
 `best-rm` means: keep the best mapping,
 if there are multiple equally good mappings,
 remove this read / contig.
 `rm` or `remove` means: remove everything that's not uniquely mapped.
 `keep` means: keep all the mappings.
 Keep in mind that with this option if a contig is mapped to multiple positions,
 a read that is mapped uniquely to it
 will end up with ambigous mapping to the reference.
 
### Example datasets

To add.

### References

Heng Li. Minimap2: pairwise alignment for nucleotide sequences;
Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

D.R. Zerbino and E. Birney. 2008. Velvet: algorithms for de novo
short read assembly using de Bruijn graphs.
Genome Research,18:821-829

