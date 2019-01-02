# IMCA

## Improving Mapping with Contigs Assembly

#### Instalation

No installation needed.
At least not for these scripts.
You might need to install some software needed to complete additional steps.
Specifically,
you need some mapper for your reads, mapper for contigs, and assembler to create contigs.
The scripts was tested using minimap2 and velvet,
but should work with output of any other software,
as long as formats are correct:
the mappings should be in .bam / .sam format,
and contigs in fasta format.

#### Workflow

1. Map reads to the reference.
2. Construct contigs from the reads.
3. Potentially, filter contigs.
4. Map (align) contigs to the reference.
5. Map reads to the mapped contigs.
6. Transfer the mapping of the reads from contigs to reference.

Steps from 1 to 5 can be done with any software of your choice.
Step 6 can be done using transferMapping.py script.
Below I present example workflow.

1. Map reads to the reference.

Using minimap2 (link):
https://github.com/lh3/minimap2

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

4. Map (align) contigs to the reference.

Using minimap2:

`minimap2 -a -x asm10 reference.fa velvet_output/contigs.fa > contigs2reference.sam 2> mapping_contigs.err`

You can also use any other mapper capable of mapping large contigs
or tool designed specifically for aligning genomes,
like MUMmer, BLAST etc.

Important note: the output of this step must be in sam or bam format.
If you use a software that generates anything else,
you need to convert it to sam / bam format before last step.
(maybe some links)
For better accuracy, CIGAR strings should be available for every alignment.
See documentation for details.

5. Map reads to the mapped contigs.

This task is simmilar to the step 1; we just use mapped and potentially filtered contigs instead of reference.
For example, using minimap2:

`minimap2 -a -x sr velvet_output/contigs.fa reads.fq > reads2contigs.sam 2> mapping2contigs.err`

But you can use Bowtie2, BWA, Tophat etc.

6. Transfer the mapping of the reads from contigs to reference.

`./transferMapping.py -r reads2reference.sam -c reads2contigs.sam -m contigs2reference.sam -d transfering_output`

This will generate directory `transfering_output`
 with files ...

#### Parameters

#### References

minimap2
https://academic.oup.com/bioinformatics/article-abstract/34/18/3094/4994778?redirectedFrom=fulltext

D.R. Zerbino and E. Birney. 2008. Velvet: algorithms for de novo
short read assembly using de Bruijn graphs.
Genome Research,18:821-829
