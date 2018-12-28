IMCA

Improving Mapping with Contigs Assembly

No installation needed.

Workflow:

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

minimap2 -a -x sr reference.fa reads.fq > reads2reference.sam 2> mapping.err

You can also use any other mapper, like bowtie2, BWA, Tophat etc.

2. Construct contigs from the reads.

Using velvet (link):

`velveth velvet_output 27 -fastq reads.fq`
`velvetg velvet_output`

This will generate directory `velvet_output` with several files,
including large Sequences, Roadmap, Graph and LastGraph, that you can remove
once `velvetg` ends running.
We are only going to need contigs.fa and potentially stats.txt for filtering.

You can use any other assembler.
Important thing is to get contigs in fasta format in this step.

3. Potentially, filter contigs.

You might choose to filter out contigs that are too short,
or have to small coverage.
I encourage to look at the obtained contigs and decide what seems best
to obtain better results.

4. Map (align) contigs to the reference.

Using minimap2:

`minimap2 -a -s asm10 reference.fa velvet_output/contigs.fa > contigs2reference.sam 2> mapping_contigs.err`

You can also use any other mapper capable of mapping large contigs
or tool designed specifically for aligning genomes,
like MUMmer, BLAST etc.
Important note: the output of this step must be in sam or bam format.
If you use a software that generates anything else,
you need to convert it to sam / bam format before last step.
(maybe some links)

5. Map reads to the mapped contigs.

This task is simmilar to the step 1, we just use mapped and potentially filtered contigs instead of reference.
For example, using minimap2:

`minimap2 -a -x sr velvet_output/contigs.fa reads.fq > reads2contigs.sam 2> mapping2contigs.err`

But you can use Bowtie2, BWA, Tophat etc.

6. Transfer the mapping of the reads from contigs to reference.

`./transferMapping.py -r reads2reference.sam -c reads2contigs.sam -m contigs2reference.sam -d transfering_output`


