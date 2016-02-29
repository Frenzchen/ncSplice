# To do

- implement single-end/paired-end option
	- filter for singletons (samtools view -f 0x08 unmmaped.bam > singletons.txt), get some kind of preparation script to prepare singletons and unmapped.fastq from unmapped.bam
	- need to check whether everything works with paired-end-option
- clean up motif scoring to trim alignments, needs to be implement in downstream analysis after seed-extension 
- Write documentation for all functions
- Re-think structure, right now
	- module Alignment: functions that manipulate data (alignments, reverse-complement, motif scoring)	
	- module Analysis: different steps of analysis
	- class ReadBam: read bam file from IO and get attributes (qname, strand, coordinates etc)
- Write overall documentation
	- versions
	- how to run
	- potential errors
- Write tests
- make gem