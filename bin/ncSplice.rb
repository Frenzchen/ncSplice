#!/usr/bin/env ruby

# 24/08/2015
# version: ruby 2.0.0
# 
# Wrapper script for circRNA detection.
# Scripts found in ../scripts, starting with circles_.
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
# 1. Handle wrong/missing/mandatory options in a better way.
# 2. Implement additional option to delete intermediate files.
# 3. Include path variable to avoid the "../scripts".


require 'optparse'
require 'open3'
require_relative "../lib/lib.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = 'Usage: ruby ruby_wrapper.rb [options]'

  opts.on('-h', '--help', 'Display help screen') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print current program versions') do
    puts "### ncSplice version ###"; puts "1.0.0"
    puts "### ruby version ###"; puts "#{RUBY_VERSION}"
    puts "### samtools version ###"; "#{system("samtools --version")}"
    puts "### bowtie2 version ###"; "#{system("bowtie2 --version")}"
    exit
  end

	options[:anchor] = 20
	
  opts.on('-i', '--input-fastq <filename>', 'input fastq file') {|i| options[:input] = i}
  opts.on('-b', '--base-name <string>', 'base name') {|b| options[:base] = b}
  opts.on('-x', '--bowtie-index <directory>', 'bowtie-index diretory and base name: <index-directory>/<bt2-index>') {|x| options[:bowtie] = x}
  opts.on('-a', '--anchor-length <integer>', Integer, 'anchor length (default 20)') {|a| options[:anchor] = a}
  opts.on('-l', '--read-length <integer>', Integer, 'read length') {|l| options[:rl] = l}  
  opts.on('-f', '--fasta-files <directory>', 'directory with fasta files, one per chromsome') {|f| options[:fasta] = f}
  opts.on('-c', '--exclude-chromosomes filename', 'text file with chromosomes to exclude') {|c| options[:exclude] = c}
end

# Remark 1
begin
  optparse.parse!(ARGV)
rescue OptionParser::InvalidOption
  $stderr.print "Error in calling scripts: #{$!}"
end


# run
##########################################################################################

# Initiate
input_file = options[:input]
base_name = options[:base]
bowtie_index = options[:bowtie]
anchor_length = options[:anchor]
read_length = options[:rl]
options[:fasta][-1] == '/' ? fasta = options[:fasta] : fasta = "#{options[:fasta]}/"
skip = options[:exclude]
logfile = File.open("logfile_#{base_name}.log", 'w')

# Preparation of anchors
begin
	Analysis.prepare_anchorpairs(input_file, anchor_length, base_name)
	Analysis.bowtie_call(bowtie_index, base_name, logfile)

	Open3.popen3("samtools view #{base_name}.bam") do |stdin, stdout, stderr, t|
		@anchor_pairs = Analysis.process_bam(stdout, fasta, skip)
	end

	Analysis.seed_extension(@anchor_pairs, anchor_length, read_length, fasta, base_name)	
	
	
rescue StandardError => err
	logfile.puts "#{Time.new}: Error in circRNAs.rb"
	logfile.puts err.message
	err.backtrace.each {|line| logfile.puts line}
	exit
end
