## Introduction


### What is ncSplice?
ncSplice is a collection of ruby scripts to identify circular, intra- and inter-chromosomal junction reads from RNASeq Illumina data. Reads indicative for these splicing events are usually discarded by standard mapping tools and therefore, are found in the unmapped fraction. For each event type, ncSplice takes the 20 most left and most right bp from each unmapped read (= anchors) and remaps them independently. Different filtering criteria are then applied to identify circular, intra- and inter-chromosomal candidate reads.

1. circular RNAs (circRNAs)
circRNAs are characterized by their atypical splicing behavior, in which the 5’-end of an exon is spliced to the 3’-end of an upstream exon (head-to-tail configuration). To be counted as a candidate read, anchor pairs have to 

	- be on the same chromosome
	- map to the same strand
	- should not be further than 100 kb apart from each other

2. Intra-chromosomal fusions
Intra-chromosomal fusion reads are defined by reads for which the anchor pair maps 

	- on the same chromosome
	- at least 1 mp apart from each other
 	- the strand orientation of anchors within the anchor pair can be different

3. Inter-chromosomal fusions
Inter-chromosomal fusion reads are defined by reads for which the anchor pair maps 

	- on different chromosomes
	- the strand orientation of anchors within the anchor pair can be different

### For what data types does ncSplice work?
Currently, the detection works for unstranded, single-end data only. Further library options (paired-end, stranded) will be implemented soon.


### Detection steps
The detection of circRNAs, intra- and inter-chromosomal fusion follows similar steps with adaptations for each of the different splice types. The general outline is:

1. Creation of anchors from unmapped reads
Take a fastq-file, which contains the unmapped reads and prepare an output fastq-file with corresponding anchors. A new qname is created, which is based on the original qname, the mate information, the read itself and the terminus A or B showing which side of the read was taken (A = left, B = right): 

```@HWI-D00108:213:C3U67ACXX:1:1106:5218:96673_1_TTTCTGTGAGCTTATGAGGCCATTCTGCACATTATCAAAATGAAATCATTATGCAGTAACCTTATACAAATCTCATAATAAAATAACCATGTATACCACT_A```

2. 1st mapping round
Anchor mapping with bowtie, output sorted according to read name.

3. Seed extension
Reads are read from the bowtie output bam-file. Each anchor pair is evaluated to see whether anchor pairs fulfill splice type conditions (see Introduction). If so, the anchors are extended. A maximum of 1 mismatch on each of the sides is allowed. If the read can be fully extended, it is written to the output file. Reads for which R1 and R2 fall on the same or different junction are removed.

4. Collision of reads on candidate loci

5. Conversion of loci into fasta format
Conversion of breakpoint into fasta-format by joining approximately 100 bp up- and downstream from the breakpoint.


6. 2nd remapping round on junction index
An index based on the detected junctions is created and all unmapped reads are remapped to find reads that span these junctions with less than 20 bp (but at least 8 bp).


7. Filtering of remapped reads and final candidate list
Uniquely remapped reads with a maximum of 2 mismatches are added to the first candidate list (if not already used).


8. A post-filtering step should be applied to remove low-coverage junctions and extreme outliers.



## ncSplice scripts, structure and system requirements

To run ncSplice, bowtie (v1 or v2) and samtools have to be added to path. They can be downloaded and installed from:

  * [bowtie1](http://bowtie-bio.sourceforge.net/bowtie1/index.shtml)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [samtools](http://samtools.sourceforge.net/)

The scripts were developed on ruby 2.0.0. The latest ruby version can be downloaded from

  * [ruby](https://www.ruby-lang.org/en/)

ncSplice will for sure not work on versions < 1.9.2. These versions do not contain the built-in method 'require_relative', which is used in ncSplice.


## How to run ncSpice

```Usage: ncSplice.rb -q <unmapped.fastq> -p <prefix> -x <index-directory>/<bt2-index> -a <anchor-length> -l <read-length> -f <chromsomes>/*.fa -c <exclude.txt> [options]```

    -h, --help                       Display help screen.
    -v, --version                    Print ncSplice version and dependencies.
    -q, --input-fastq <filename>     Unmapped reads in fastq format.
        --sequencing-type <string>   Sequencing type, 'se' for single-end and 'pe' for paired-end sequencing, default is 'se'
        --singletons <string>        Single mapped reads, only required if sequencing-type set to 'pe'
    -p, --prefix <string>            Prefix for all files.
    -x, --bowtie-index <directory>   Bowtie-index diretory and base name: <index-directory>/<bt2-index>.
    -a, --anchor-length <integer>    Length of the read anchor for remapping, default is 20 bp, shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in candidates.
    -l, --read-length <integer>      Length of the sequencing read.
    -f, --fasta-files <directory>    Directory with chromosome fasta files, one fasta-file per chromsome.
    -s, --skip-chr <filename>        Text file with chromosomes to exclude, such as the mitochondrial chromosome (recommanded), chromosomes need to be listed in a separate text-file with one chromosome per line.


## To do

1. Implement single-end/paired-end option
  * filter for singletons (samtools view -f 0x08 unmmaped.bam > singletons.txt), get some kind of preparation script to prepare singletons and unmapped.fastq from unmapped.bam
  * need to check whether everything works with paired-end-option
2. Clean up motif scoring to trim alignments, needs to be implement in downstream analysis after seed-extension 
3. Write documentation for all functions
4. Re-think structure, right now
  * module Alignment: functions that manipulate data (alignments, reverse-complement, motif scoring)
  * module Analysis: different steps of analysis
  * class ReadBam: read bam file from IO and get attributes (qname, strand, coordinates etc)
5. Write overall documentation
  * versions
  * how to run
  * potential errors
6. Write tests
7. Make gem