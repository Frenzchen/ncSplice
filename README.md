# Project status

- subscripts into functions, stored in lib.rb

# To do

- implement single-end/paired-end option
- clean up motif scoring to trim alignments, needs to be implement in downstream analysis after seed-extension 
- Work on log-file
- Clean up option-passing
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