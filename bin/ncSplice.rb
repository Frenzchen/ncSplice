#!/usr/bin/env ruby

# 13/03/2016
# version: ruby 2.0.0
# 
# Wrapper script for circRNA detection.
#
#
# Options (update!):
# u - Fastq-file with unmapped reads from tophat in bam-format.
# p - Prefix used for all output files.
# x - Path to bowtie2-index directory and base name <index_directory>/<bt2-index>.
# a - Length of anchors, 20 bp is recommended.
# l - Length of read.
# f - Path to directory with one fasta files for each chromosome.
# s - File with chromosomes that should excluded.
#			One chromosome per line.
# singletones - File with paired reads for which only one mate mapped.
# sequencing-type - SE or PE for single-end or paired-end.
#
# Remarks:
# Implement additional option to delete intermediate files.


require 'optparse'
require 'open3'
require_relative "../lib/alignments.rb"
require_relative "../lib/analysis.rb"
require_relative "../lib/bamClass.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@isb-sib.ch)\n
  Usage: ncSplice.rb -q <unmapped.fastq> -p <prefix> -x <index-directory>/<bt2-index> -a <anchor-length> -l <read-length> -f <chromsomes>/*.fa -c <exclude.txt> [options]"

  opts.on('-h', '--help', 'Display help screen.') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print ncSplice version and dependencies.') do
		puts '# ncSplice'
		puts 'ncSplice v0.1.0'
		puts "\n# Dependencies"
		puts ['ruby >= 2.0.0', 'samtools >= 1.0.0', 'bowtie2 >= 2.1.0'].join("\n")
    exit
  end

	options[:phred] = 25
	options[:anchor] = 20
	options[:sequencing] = 'se'
	
	opts.on('-u', '--unmapped <filename>', String, 'Bam file with unmapped reads') {|u| options[:unmapped] = u}
	opts.on('-q', '--quality <integer>', Integer, 'Minimal phred quality unmapped reads need to have for further analysis.') {|q| options[:phred] = q}
  opts.on('--sequencing-type <string>', String, 'Sequencing type, \'se\' for single-end and \'pe\' for paired-end sequencing, default is \'se\'') {|seq| options[:sequencing] = seq}
  opts.on('--singletons <string>', String, 'Single mapped reads, only required if sequencing-type set to \'pe\'') {|seq| options[:singletons] = seq}
  opts.on('-p', '--prefix <string>', String, 'Prefix for all files.') {|b| options[:prefix] = b}
  opts.on('-x', '--bowtie-index <directory>', String, 'Bowtie-index diretory and base name: <index-directory>/<bt2-index>.') {|x| options[:bowtie] = x}
  opts.on('-a', '--anchor-length <integer>', Integer, 'Length of the read anchor for remapping, default is 20 bp, shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in candidates.') {|a| options[:anchor] = a}
  opts.on('-l', '--read-length <integer>', Integer, 'Length of the sequencing read.') {|l| options[:readlength] = l}  
  opts.on('-f', '--fasta-files <directory>', 'Directory with chromosome fasta files, one fasta-file per chromsome.') {|f| options[:fasta] = f}
  opts.on('-s', '--skip-chr <filename>', String, 'Text file with chromosomes to exclude, such as the mitochondrial chromosome (recommanded), chromosomes need to be listed in a separate text-file with one chromosome per line.') {|s| options[:skip] = s}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length <= 6

# local
unampped_bam = options[:unmapped]
phred_quality = options[:phred]
anchor_file = options[:input]
prefix = options[:prefix]
bowtie_index = options[:bowtie]
anchor_length = options[:anchor]
read_length = options[:readlength]
options[:fasta][-1] == '/' ? fasta = options[:fasta] : fasta = "#{options[:fasta]}/"
skip = options[:skip]
singletons = options[:singletons]

# global
$sequencing_type = options[:sequencing]
$logfile = File.open("#{prefix}_logfile.log", 'w')


# run
##########################################################################################

begin
	Open3.popen3("samtools view #{unmapped_bam}") do |stdin, stdout, stderr, t|
		Analysis.unmapped2fastq(stdout, "#{prefix}_unmapped.fastq", phred_quality)
	end
	
	if $sequencing_type == 'pe'
		$singletons = Analysis.read_singletons(singletons, read_length) 
	end
	
	Analysis.prepare_anchorpairs("#{prefix}_unmapped.fastq", anchor_length, "#{prefix}_anchors.fastq")
	Analysis.bowtie_map(bowtie_index, "#{prefix}_anchors.fastq", "#{prefix}.bam")

	anchor_pairs = Open3.popen3("samtools view #{prefix}.bam") do |stdin, stdout, stderr, t|
		Analysis.process_bam(stdout, fasta, skip)
	end

	Analysis.seed_extension(anchor_pairs, anchor_length, read_length, fasta, "#{prefix}_candidateReads.txt")	
	Analysis.collaps_qnames("#{prefix}_candidateReads.txt", "#{prefix}_candidates.txt")
	Analysis.candidates2fa("#{prefix}_candidates.txt", fasta, read_length, "#{prefix}_faIndex.fa")
	Analysis.bowtie_build("#{prefix}_faIndex.fa", "#{prefix}")
	
	Open3.popen3("mkdir #{prefix}_index") if !Dir.exists?("#{prefix}_index")
	Open3.popen3("mv #{prefix}_faIndex.fa #{prefix}_index/; mv *bt2 #{prefix}_index/")
	
	Analysis.bowtie_map("#{prefix}_index/#{prefix}", anchor_file, "#{prefix}_remapping.bam")

	Open3.popen3("samtools view #{prefix}_remapping.bam") do |stdin, stdout, stderr, t|
		Analysis.remapped_reads(stdout, "#{prefix}_remapped.txt", read_length)
	end
	
	Analysis.final_candidates("#{prefix}_candidates.txt",  "#{prefix}_remapped.txt", "#{prefix}_final.txt")
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
