#!/usr/bin/env ruby

# 28/02/2016
# version: ruby 2.0.0
# 
# Wrapper script for circRNA detection.
#
#
# Options:
# i - Fastq-file with unmapped reads.
# b - Base name with which all output-files will start.
# x - Path to bowtie2-index directory and base name <index_directory>/<bt2-index>.
# a - Length of anchors, 20 bp is recommended.
# l - Length of read.
# f - Path to directory with one fasta files for each chromosome.
# c - File with chromosomes that should excluded.
#			One chromosome per line.
#
# Remarks:
# 2. Implement additional option to delete intermediate files.


require 'optparse'
require 'open3'
require_relative "../lib/lib.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@unil.ch)\n
  Usage: ncSplice.rb -q <unmapped.fastq. -p <prefix> -x <index-directory>/<bt2-index> -a anchor-length -l read-length -f <chromsomes>/*.fa -c exclude.txt"

  opts.on('-h', '--help', 'Display help screen') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print dependencies') do
	  puts "# ncSplice"
    puts "ncSplice v0.1.0"
		puts "\n# Dependencies"
    puts ["ruby >= 2.0.0", "samtools >= 1.0.0", "bowtie2 >= 2.1.0"].join("\n")
    exit
  end

	options[:anchor] = 20
	options[:sequencing] = 'se'
	
  opts.on('-q', '--input-fastq <filename>', String, 'Unmapped reads in fastq format.') {|q| options[:input] = q}
  #opts.on('--sequencing-type <string>', String, 'Sequencing type, se for single-end and pe for paired-end sequencing') {|seq| options[:sequencing] = seq}
  opts.on('-p', '--prefix <string>', String, 'Prefix for all files.') {|b| options[:prefix] = b}
  opts.on('-x', '--bowtie-index <directory>', String, 'Bowtie-index diretory and base name: <index-directory>/<bt2-index>.') {|x| options[:bowtie] = x}
  opts.on('-a', '--anchor-length <integer>', Integer, 'Length of the read anchor for remapping, default is 20 bp, shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in candidates.') {|a| options[:anchor] = a}
  opts.on('-l', '--read-length <integer>', Integer, 'Length of the sequencing read.') {|l| options[:readlength] = l}  
  opts.on('-f', '--fasta-files <directory>', 'Directory with chromosome fasta files, one fasta-file per chromsome.') {|f| options[:fasta] = f}
  opts.on('-s', '--skip-chr <filename>', String, 'Text file with chromosomes to exclude, such as the mitochondrial chromsome MT (recommanded), chromsomes need to be listed in separate text-file with one chromosome per line.') {|s| options[:skip] = s}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length <= 6

# local
input_file = options[:input]
prefix = options[:prefix]
bowtie_index = options[:bowtie]
anchor_length = options[:anchor]
read_length = options[:readlength]
options[:fasta][-1] == '/' ? fasta = options[:fasta] : fasta = "#{options[:fasta]}/"
skip = options[:skip]

# global
$sequencing_type = options[:sequencing]
$logfile = File.open("#{prefix}_logfile.log", 'w')


# run
##########################################################################################

begin
	$singletons = Analysis.read_singletons(singletons) if $sequencing_type == 'pe'

	Analysis.prepare_anchorpairs(input_file, anchor_length, "#{prefix}_anchors.fastq")
	Analysis.bowtie_map(bowtie_index, "#{prefix}_anchors.fastq", "#{prefix}.bam")

	anchor_pairs = Open3.popen3("samtools view #{prefix}.bam") do |stdin, stdout, stderr, t|
		Analysis.process_bam(stdout, fasta, skip)
	end

	Analysis.seed_extension(anchor_pairs, anchor_length, read_length, fasta, "#{prefix}_candidateReads.txt")	
	Analysis.collaps_qnames("#{prefix}_candidateReads.txt", "#{prefix}_candidates.txt")
	Analysis.candidates2fa("#{prefix}_candidates.txt", fasta, read_length, "#{prefix}_faIndex.fa")
	Analysis.bowtie_build("#{prefix}_faIndex.fa")
	
	Open3.popen3("mkdir #{prefix}_index") if !Dir.exists?("#{prefix}_index")
	Open3.popen3("mv #{prefix}_faIndex.fa #{prefix}_index/; mv *bt2 #{prefix}_index/")
	
	Analysis.bowtie_map("index_circles/candidates", input_file, "#{prefix}_remapping.bam")

	Open3.popen3("samtools view #{prefix}_remapping.bam") do |stdin, stdout, stderr, t|
		Analysis.remapped_reads(stdout, "#{prefix}_remappedCandidates.txt", read_length)
	end
	
	Analysis.final_candidates("#{prefix}_candidates.txt",  "#{prefix}_remappedCandidates.txt", "#{prefix}_final.txt")
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
