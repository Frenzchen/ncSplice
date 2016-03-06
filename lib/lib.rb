#!/usr/bin/env ruby

module Alignment

	extend self
	
	# Trim end of alignment if mismatch is in the last 2 bp, likely to be a mapping error.
	#
	#
	# array	- Array of Integers from mapping with breakpoint_downstream or
	#					breakpoint_upstream methods.
	# mm		- Number of allowed mismatches.
	#
	# Return integer of trimmed alignment array.
	def trim(array, mm)
		array = array[0..-3] if array[-2..-1] == [0,1] && mm == 1
		array = array[0..-2] if array[-1] == 0 && mm == 1
		array.length
	end


	# Extend 3' anchor to upstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def upstream(seq, reference, mm)
		a = []
		mismatch = 0

		(seq.length - 1).downto(0).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Extend 5' anchor to downstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def downstream(seq, reference, mm)
		a = []
		mismatch = 0
	
		0.upto(seq.length - 1).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Create reverse complement of DNA.
	#
	#
	# dna - DNA sequence
	#
	# Return reverse complement.
	def reverse_complement(dna)
		complement = []
		dna.each_char do |s|
			if s == 'A'
				complement << 'T'
			elsif s == 'T'
				complement << 'A'
			elsif s == 'G'
				complement << 'C'
			elsif s == 'C'
				complement << 'G'
			elsif s == 'N'
				complement << 'N'
			end	
		end
		complement.join('').reverse
	end
	
	
	
	def quality_ok?(quality, phred)
  	quality.each_char.all? {|x| (x.ord - 33) >= phred}
	end



	# !!! Function still under development !!!
	#
	# Compare motif via 2bp sliding window to know U2 and U12 splice sites.
	# Each motif pair (upstream + downstream) is scored.
	# score += 0.5 if motif is present in acceptor or donor motis.
	# Motif pair with maximum score is defined and the index of the motof returned
	#
	# Returns integer corrsponding to index of highest scored motif
	def score_motifs(upstream_dna, downstream_dna, overhang)
		splicesites = {:acc => ['GT', 'AT'], :don => ['AG', 'AC']}
		l = overhang - 1
		fwd_scores = {:acc => Array.new(l, 0), :don => Array.new(l, 0)}
		rev_scores = {:acc => Array.new(l, 0), :don => Array.new(l, 0)}
		up_fwd = upstream_dna[0..overhang-1]
		down_fwd = downstream_dna[-overhang..-1]
		up_rev = reverse_complement(upstream_dna[0..overhang-1])
		down_rev = reverse_complement(downstream_dna[-overhang..-1])

		# scores for sense strand
		0.upto(overhang-1) do |i|
			fwd_scores[:don][i] += 0.5 if splicesites[:don].include?(up_fwd[i..i+1])	
			fwd_scores[:acc][i] += 0.5 if splicesites[:acc].include?(down_fwd[i..i+1])
		end
		
		# scores for reverse strand
		0.upto(overhang-1) do |i|
			rev_scores[:acc][i] += 0.5 if splicesites[:acc].include?(up_rev[i..i+1])	
			rev_scores[:don][i] += 0.5 if splicesites[:don].include?(down_rev[i..i+1])
		end
				
		# find maximal motif score	
		summed_scores_fwd = (fwd_scores.values.transpose.map {|x| x.reduce(:+)})
		summed_scores_rev = (rev_scores.values.transpose.map {|x| x.reduce(:+)})
		max_fwd = summed_scores_fwd.max 
		max_rev = summed_scores_rev.max
		
		if max_fwd > 0 || max_rev > 0 

			# fwd strand with higer score
			if max_fwd > max_rev 
				i = summed_scores_fwd.index {|x| x == max_fwd}
				o = overhang - (2*i)-1
				motif = [up_fwd[i..i+1], down_fwd[i+2..-1]]
				{:index => i, :score => max_fwd, :strand => 1, :motif => motif.join('|'), \
				 :overhang => o}
			
			# rev strand with higher score
			elsif max_fwd < max_rev
				o = overhang - 2*(i+1)
				i = summed_scores_rev.index {|x| x == max_rev}	
				motif = [up_rev[i..i+1], down_rev[i+2..-1]]
				{:index => i, :score => max_rev, :strand => -1, :motif => motif.join('|'), \
				 :overhang => o}
			
			# equal scores on fwd and rev strand
			else
		 		i_fwd = summed_scores_fwd.index {|x| x == max_fwd}
		 		i_rev = summed_scores_rev.index {|x| x == max_rev}
	
				# same/different index position
		 		if i_fwd == i_rev
		 			i = i_fwd
		 			o = overhang - 2*(i+1)
		 			motif = [up_fwd[i..i+1], down_fwd[i+2..-1]]
		 			{:index => i, :score => max_fwd, :strand => 0, :motif => motif.join('|'), \
		 			 :overhang => o}
		 		else
		 			{:index => '-', :score => max_fwd, :strand => 0, :motif => '-', \
		 			 :overhang => overhang}
		 		end
		 	end
		# unknown motif	
		else
		 	{:index => '-', :score => '-', :strand => 0, :motif => '-', :overhang => overhang}
		end
	end


	# Compare number of mismatches reported by MD:Z tag to allowed number of mismatches.
	#
	#
	# mdz - String of MD:Z tag from bowtie2 alignment.
	# mm  - Integer of max. number of allowed mismatches.
	#
	# Return boolean.
	def mismatches?(mdz, mm)
		mdz = mdz.split(':').last
		counter = 0
		%w[A T G C].each {|x| counter += mdz.count(x) if !mdz.nil?}
		counter > mm
	end
	
	
	# Report genompic distance over which read mapped.
	#
	#
	# cigar				- Cigar string from bowtie2 alignment.
	# read_length	- Length of sequecning reas.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def genomic_mappinglength(cigar, read_length)
		if cigar == '{read_length}M'
			read_length
		elsif cigar.include?('N') && %w(I D S H P).all? {|x| !cigar.include?(x)}
			cigar = cigar.split(/N|M/).collect {|x| x.to_i}
			cigar.inject {|sum, x| sum + x}
		else
			false
		end
	end
	
	
	# Check whether spliced read has mapped mate in the right position and orientation.
	#
	#
	# mate1	- Mapping information for mate 1.
	# mate2	- Mapping information for mate 2.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def paired?(mate1, mate2)
		mate1_chr, mate1_start, mate1_stop, mate1_strand = mate1
		mate2_chr, mate2_start, mate2_stop, mate2_strand = mate2
		conditions = []
		
		conditions.push(mate1_chr == mate2_chr)
		conditions.push(mate1_start <= mate2_start && mate1_stop >= mate2_stop)
		conditions.push(mate1_strand != mate2_strand)
		conditions.all? { |con| con == true }
	end
end

##########################################################################################

class ReadBam

	# Holds methods for backsplice, intra- and inter-chromosomal junction detection. 
	# Usable for unstranded, paired-end sequencing data.
	# Takes sam/bam line to create object.
	#
	#
	# Create new SplicedRead object:
	#
	# id      - Read id of the read anchor.
	# strand  - Strand on which the anchor maps. 
	# chr     - Chromosome on which the anchor maps.
	# start   - Start position of anchor mapping.
	# cigar   - Cigar string of anchor mapping.
	attr_accessor :id, :strand, :chr, :start, :cigar
	
	def initialize(line)
		@id, @strand, @chr = line[0..2]
		@start = line[3].to_i
		@cigar = line.find {|l| l.match(/MD:Z:[[:digit:]]+/) }
	end
	
	def start
		@start.to_i
	end
	
	def strand
		@strand.to_i & 0x10 > 0 ? @strand = -1 : @strand = 1
	end
	
	
	# Examines whether read pairs fulfills circular RNA conditions.
	# Those are:
	# 1. Anchors have to be on the same strand.
	# 2. Anchors need to be on the same chromosome.
	# 3. Anchors need to match completely.
	# 4. Chromosome should not be on the "to exclude"-list.
	# 5. Anchors should map a max. of 100 kb apart from each other.
	# 6. Anchors have to be in head-to-tail configuration.
	#
	# anchor		- Other half of anchor pair, given as object created by SplicedRead class.
	# distance	- Integer of the max. distance anchors should be apart form each other.
	# exclude		- Array with chromosome names to ignore.
	#
	# Returns boolean.
	def circular?(anchor, distance, exclude)
		conditions = []
		
		conditions << (strand == anchor.strand) # 1.
		conditions << (chr == anchor.chr) # 2.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'} # 3.
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)}	# 4.
		conditions << ((start - anchor.start).abs <= distance) # 5.
		if strand == 1 && anchor.strand == 1 # 6.
			conditions << (start > anchor.start)
		elsif strand == -1 && anchor.strand == -1
			conditions << (start < anchor.start)
		end
	
		conditions.all? {|a| a == true}
	end
	
	
	# Examines whether read pairs fulfills intra-chromosomal transcript conditions.
	# Those are:
	# 1. Anchors need to be on the same chromosome.
	# 2. Anchors need to match completely.
	# 3. Chromosome should not be on the "to exclude"-list.
	# 4. Anchors should map a min. of 1 mb apart from each other.
	#
	# anchor		- Other half of anchor pair given as object created by SplicedRead class.
	# distance	- Integer of the max. distance anchors should be apart form each other.
	# exclude		- Array with chromosome names to ignore
	#
	# Returns boolean.
	def intraChimeric?(anchor, distance, exclude)
		conditions = []
		
		conditions << (chr == anchor.chr) # 1.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'} # 2.
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)}	# 3.
		conditions << ((start - anchor.start).abs >= distance) # 4.
	
		conditions.all? {|a| a == true }
	end


	# Examines whether read pairs fulfills intra-chromosomal transcript conditions.
	# Those are:
	# 1. Anchors need to be on different chromosomes.
	# 2. Anchors need to match completely.
	# 3. Chromosome should not be on the "to exclude"-list.
	#
	#	anchor	- Other half of anchor pair given as object created by SplicedRead class.
	# exclude	- Array with chromosome names to ignore.
	#
	# Returns boolean.
	def interChimeric?(anchor, exclude)
		conditions = []
		
		conditions << (chr != anchor.chr) # 1.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'}	# 2.	
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)} # 3.

		conditions.all? {|a| a == true}
	end
end


##########################################################################################

module Analysis

	extend self
	
	# Catch exit status of system process and report it to logfile.
	# Exit program if subprocess was not finished successfully.
	#
	#
	# t       - Open3 object for thread.
	# stderr  - Open3 object for STDERR.
	# name    - Name for logfile.
	# logfile - Name of logfile to write to.
	#
	# Return string that is written to logfile. 
	def system_exitcode(t, stderr, name)	
		if t.value.success?
			$logfile.puts "#{Time.new.strftime("%c")}: #{name} succeded."
			if stderr.any?
				$logfile.puts "#{name} output:"
				stderr.readlines.each {|line| $logfile.puts line}
			end
		else
			$logfile.puts "#{Time.new.strftime("%c")}: Error in #{name}:"
			stderr.readlines.each {|line| $logfile.puts line}
			exit
		end
	end
 
 
  # Convert unmapped.bam to fastq-format
 	# Do directly from accepted_hits.bam or from separate file??
 	#
 	#
 	# 
 	# 
 	#
 	# 
 	def unmapped2fastq(input_file, output_file, phred_qualit)
 		File.open(output_file, 'w') do |output|
			input_file.each do |line|
  			line = line.strip.split(/\s+/)
  
  			flag = line[1].to_i
  			flag & 0x40 > 0 ? mate = '1' : mate = '2'
  	
  			qname, sequence, quality = line[0], line[9], line[10] 
  			output.puts "@#{qname}/#{mate}", sequence, '+', quality if quality_ok?(quality, phred_quality)
  		end
  	end
  	$logfile.puts "#{Time.new.strftime("%c")}: Converted unmapped.bam into fastq-format."	
	end

 	
 	# Read singletons.
 	# Do directly from accepted_hits.bam or from separate file??
 	#
 	#
 	# singletons	- Text-file with singletons (keep as bam-file?)
 	# read_length	- Length of sequencing read.
 	#
 	# Return hash with all singletons, key = qname, values = mapping information.
	def read_singletons(singletons, read_length)
		single_reads = {}
		
		File.open(singletons, 'r').readlines.each do |line|
  
  		line = line.strip.split(/\s+/)
  		qname, flag, chr, start = line[0..3]  	
  		flag.to_i & 0x10 > 0 ? strand = -1 : strand = 1
  		cigar = line[5]
			distance = genomic_mappinglength(cigar, read_length)
			
			if distance != false
				strand == 1 ? stop = start + distance : stop = start - distance
				single_reads[qname] = [chr, start, stop, strand]
			end
		end
		single_reads
	end
	
	
	# Prepare fastq-file with anchors from unmapped.fastq.
	#
	#
	# input_file    - unmapped.fastq.
	# anchor_length - Length of anchor, default is 20 bp.
	# output_file   - Name of output file.
	#
	# Return fastq-file with anchor pairs.
	def prepare_anchorpairs(input_file, anchor_length, output_file)	
		name, mate, seq, quality = nil
		counter = -1

		File.open(output_file, 'w') do |output| 
			File.open(input_file, 'r').readlines.each do |line|
				counter += 1
				line = line.strip
			
				if counter % 4 == 0 
					name, mate = line.split('/')

				elsif counter % 4 == 1
					seq = line
				
				elsif counter % 4 == 3
					quality = line
			
					name_A = "#{name}_#{mate}_#{seq}_A"
					name_B = "#{name}_#{mate}_#{seq}_B"
			
					seq_A = seq[0..anchor_length - 1]
					seq_B = seq[-anchor_length..-1]
	
					quality_A = quality[0..anchor_length - 1]
					quality_B = quality[-anchor_length..-1]
			
					output.puts [name_A, seq_A, '+', quality_A, name_B, seq_B, '+', quality_B].join("\n")
				
					name, mate, seq, quality = nil
					counter = -1
				end 
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Anchor preparation succeded."	
	end


	# Mapping of anchors via bowtie2 system call.
	#
	#
	# bowtie_index  - Path to bowtie index and index base name
	# fastq_file    - Fastq-file with anchors
	# output_file   - Name of output file
	# logfile       - Name of logfile to write to.
	#
	# Return mapped anchors in bam-format.
	def bowtie_map(bowtie_index, fastq_file, output_file)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -q -U #{fastq_file} | samtools view -bS - > #{output_file}")
		system_exitcode(t, stderr, 'Bowtie2')
	end


	# Read bam-file from IO and process lines pair-wise to filter circRNA candidates.
	#
	#
	# input_file - Input file from IO
	# fasta      - Path to chromsomal fasta-files.
	# skip       - Text-file with chromsomes to skip, one/line.
	#
	# Return hash with valid anchor pairs.
	def process_bam(input_file, fasta, skip, distance=100000)

		# general settings
		exclude = []
		File.open(skip, 'r').readlines.each {|line| exclude << line.strip}
		firstline = TRUE 
		anchor_left = nil
		anchor_right = nil
		input_hash = {}

		# Initiate chromosome hash
		Dir.foreach(fasta) do |item|
			chr = item.sub('.fa', '')
			next if item == '.' || item == '..' || exclude.include?(chr)  
			input_hash[chr] = []
		end

		# read bam file
		input_file.each do |line|
			line = line.strip.split(/\s+/)
		
			if firstline 
				anchor_left = ReadBam.new(line)
				firstline = FALSE
			else
				anchor_right = ReadBam.new(line)
				chr = anchor_right.chr

				if anchor_left.circular?(anchor_right, distance, exclude)
		
					# store coordinate sorted
					if anchor_left.strand == 1
						input_hash[chr] << [anchor_right, anchor_left] 
					else
						input_hash[chr] << [anchor_left, anchor_right]
					end
				end
	
				anchor_left, anchor_right = nil
				firstline = TRUE
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Found anchor pairs."		
		input_hash
	end


	# Seed extension for anchor pairs. They are kept as candidate if extension is succesful.
	#
	#
	# input_hash    - Input hash with anchor pairs.
	# anchor_length - Length of anchors
	# read_length   - Length of sequencing read.
	# fasta         - Path to chromsomal fasta-files.
	#	output_file   - Name of output_file.
	#
	# Return tab-delimited file with initial candidates.
	def seed_extension(input_hash, anchor_length, read_length, fasta, output_file, mm = 1, max_overhang = read_length + 8)

		output_hash = {}
	
		input_hash.each do |chr, anchorpairs|
			chr = chr.to_s
	
			# Load reference
			fasta_file = File.open("#{fasta}#{chr}.fa", 'r')
			header = fasta_file.gets.strip
			dna = fasta_file.read.gsub(/\n/, '')
	
			# Loop through hash to extend seeds for each pair
			anchorpairs.each do |pair|
				upstream, downstream = pair
				qname, mate, read = upstream.id.split('_')[0..2]

				strand, chr, upstream_start, downstream_start = upstream.strand, upstream.chr, upstream.start - 1, downstream.start - 1

				read = Alignment.reverse_complement(read) if strand == -1
				upstream_dna = dna[upstream_start - read_length + anchor_length..upstream_start + anchor_length - 1].upcase
				downstream_dna = dna[downstream_start..downstream_start + read_length - 1].upcase				
				upstream_alignmentlength = Alignment.upstream(read, upstream_dna, mm)
				downstream_alignmentlength = Alignment.downstream(read, downstream_dna, mm)
				total_alignmentlength = upstream_alignmentlength + downstream_alignmentlength

				if total_alignmentlength >= read_length && total_alignmentlength <= max_overhang
					upstream_breakpoint = upstream_start - upstream_alignmentlength + anchor_length	
					downstream_breakpoint = downstream_start + downstream_alignmentlength - 1
					overhang = total_alignmentlength - read_length

					# splice site analysis
					up = upstream_dna[-upstream_alignmentlength..-1]
					down = downstream_dna[0..downstream_alignmentlength-1]	
					motif_summary = Alignment.score_motifs(up, down, overhang)

					if motif_summary[:index].kind_of?(Integer)
						i = motif_summary[:index]
						upstream_breakpoint = upstream_breakpoint + i + 2
						downstream_breakpoint = downstream_breakpoint - overhang + i + 2
					end
	
					qname = qname.to_sym
					summary = [chr, upstream_breakpoint, downstream_breakpoint, strand, total_alignmentlength, motif_summary[:score], motif_summary[:strand], motif_summary[:motif], read_length + motif_summary[:overhang], mate] 

					# Candidates for which both, R1 and R2, are present are deleted
					# One read can neither fall on two different non-canonical nor the same junction
					if !output_hash.has_key?(qname)
						output_hash[qname] = summary
					else
						output_hash.delete(qname)
					end
				end
			end
		end
		
		# select reads with paired mate
		if $sequencing_type == 'pe'
			output_hash[qname].select! {|qname, values| $singletons.has_key?(qname) && paired?(summary, $singletons[qname])}
		end

		File.open(output_file, 'w') do |output|
			output_hash.each do |qname, v| 
				output.puts ["#{qname.to_s}/#{v[-1]}", v[0..-2]].join("\t") if v[2] - v[1] >= read_length
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Seed extension succeded."
	end
	
	
	# Collaps candidates reads onto common loci.
	#
	#
	# input_file	- Input file.
	# output_file	- Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def collaps_qnames(input_file, output_file)
	
		loci = {}
	
		# Read candidate loci and count reads/locus
		File.open(input_file, 'r').readlines.each do |line|
			line = line.strip.split("\t")
			qname = line[0]
			base = qname.gsub(/\/[1,2]/, '')
			pos = line[1..3].join(':')
			alignment_length = line[5]
	
			if !loci.has_key?(pos)
				loci[pos] = {:count => 1, :qnames => [qname], :l => alignment_length}
			else 
				loci[pos][:qnames] << qname
				loci[pos][:count] += 1
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			loci.each do |pos, v| 
				output.puts [pos.split(':'), v[:count], v[:l], v[:qnames].join(';')].join("\t") if v[:count] > 0
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Collapsed anchor pairs to single loci."
	end
	
	
	# Transform candidate loci into fasta-format.
	#
	#
	# input_file  - Input file.
	# fasta       - Path to chromsomal fasta-files.
	# read_length - Length of sequencing read.
	# output_file	- Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def candidates2fa(input_file, fasta, read_length, output_file, exoncov=8)
		chromosomes = {}
		positions = []
		
		# Input into hash sorted by chromosomes
		File.open(input_file, 'r').readlines.each do |line|
			line = line.strip.split("\t")[0..-2]
			pos = line[0..2].join(':')
			chr = line[0]
	
			if !chromosomes.has_key?(chr)
				chromosomes[chr] = [line]
		
			# 2nd elsif to exclude reads that map on same junction but opposite ends		
			elsif chromosomes.has_key?(chr) && !positions.include?(pos)
				chromosomes[chr].push(line)
				positions << pos
			end
		end

		# Output
		output = File.open(output_file, 'w') do |output|
			chromosomes.each do |chr, values|
				fasta = File.open("#{fasta}#{chr}.fa", 'r')
				header = fasta.gets.strip
				dna = fasta.read.gsub(/\n/, '')
		
				values.each do |v|
					bp_a, bp_b = v[1..2].collect {|x| x.to_i}
					overlap = v[-1].to_i - read_length
					l = read_length - exoncov 
			
					upstream = dna[bp_a..bp_a + overlap + l - 1].upcase	
					downstream = dna[bp_b - l - overlap + 1..bp_b - overlap].upcase
			
					output.puts [">#{v[0..2].join(':')}", downstream + upstream].join("\n")
				end
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Wrote loci to fasta-file."
	end
	
	
	# Built bowtie-index from candidate loci.
	#
	#
	# input_file - Input file.
	# logfile    - Name of logfile to write to.
	#
	# Return tab-delimited file with candidate loci.
	def bowtie_build(input_file, prefix)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2-build -q -f #{input_file} #{prefix}")
	system_exitcode(t, stderr, 'Building bowtie2 index')
	end
	
	
	# Filter remapped reads, keep those that haven't been used yet and fullfill criteria
	#
	#
	# input_file  - Input file.
	# output_file - Name of output_file.
	# read_length - Length of sequencing read.
	# mm          - Number of mismatchs, default 2.	
	#
	# Return tab-delimited file with candidate loci.
	def remapped_reads(input_file, output_file, read_length, mm=2)
		remapped = {}
		
		# Filter remapped reads
		input_file.each do |line|
			mdz = line.match(/MD:Z:[[:digit:]]+/).to_s
			line = line.strip.split(/\s+/)
			qname, mate = line[0].split('/')
			chr, start, stop = line[2].split(':')
			cigar = line[5]
	
			if !remapped.has_key?(qname) && !Alignment.mismatches?(mdz, mm) && cigar == "#{read_length}M"
				remapped[qname] = [chr, start, stop, mate]
			else	
				remapped.delete(qname)
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			remapped.each {|k, v| output.puts ["#{k}/#{v[-1]}", v[0..2]].join("\t")}
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Found remapped reads."
	end
	
	
	# Compare first mapping and remapping to get final candidates.
	#
	#
	# before      - Input file.
	# after       - Name of logfile to write to.
	# output_file - Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def final_candidates(before, after, output_file)
		circles = {}
		all_ids = {}

		# Read circular candidates into hash
		File.open(before, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			pos = line[0..2].join(':')
			read_count = line[3].to_i
			qname = line[-1].split(';')
	
			# Create qname index to make search faster
			# Remark 2
			qname.each do |q|
				k1, k2 = q.split(':')[3..4]
		
				all_ids[k1] = {} if !all_ids.has_key?(k1)
		
				if !all_ids[k1].has_key?(k2)
					all_ids[k1][k2] = [q]
				else
					all_ids[k1][k2] << q
				end
			end

			circles[pos] = {:counts => read_count, :qnames => qname}
		end

		# Read remapped readpairs and compare them to initial candidates
		File.open(after, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			qname = line[0]
			pos = line[1..3].join(':')
			k1, k2 = qname.split(':')[3..4]

			read_unused = (!all_ids.has_key?(k1) || !all_ids[k1].has_key?(k2) || !all_ids[k1][k2].include?(qname)) 
			
			if $sequencing_type == 'pe'
				qname = qname.split('/').first
				next if !$singletons.has_key?(qname) || !paired?(summary, $singletons[qname])
			end
			
			# Add read if read is not already used (condition 2)
			if circles.has_key?(pos) && read_unused 
				circles[pos][:counts] += 1
				circles[pos][:qnames] << qname
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			output.puts %w(chr pos_a pos_b readCounts qnames).join("\t")
	
			circles.each do |pos, v| 
				output.puts [pos.split(':'), v[:counts], v[:qnames].join(';')].join("\t")
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Final candidate list finished."
	end
end