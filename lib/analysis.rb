#!/usr/bin/env ruby

require_relative "alignments.rb"

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
 	#
 	#
 	# input_file		- Unmapped reads from tophat in bam format.
 	# output_file		-	Name of output file (in fastq).
 	# phred_quality	- Minimal phredquality for read to be kept.
 	#
 	# Unmapped reads in fastq format.
 	def bam2fastq(input_file, output_file, phred_quality)
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
	
	# Mapping of anchors via bowtie2 system call.
	#
	#
	# bowtie_index  - Path to bowtie index and index base name
	# fastq_file    - Fastq-file with anchors
	# output_file   - Name of output file
	# logfile       - Name of logfile to write to.
	#
	# Return mapped anchors in bam-format.
	def bowtie_map2(bowtie_index, fastq_file, output_file)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -q -U #{fastq_file} --no-unal | samtools view -bS - > #{output_file}")
		system_exitcode(t, stderr, 'Bowtie2')
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
end