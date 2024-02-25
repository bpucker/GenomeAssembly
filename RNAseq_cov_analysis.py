### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

#coverage file construction taken from: Pucker & Brockington, 2018: https://doi.org/10.1186/s12864-018-5360-z

__usage__ = """
	python3 RNAseq_cov_analysis.py
	--bam <BAM_INPUT_FILE> | --cov <COVERAGE_INPUT_FILE>
	--gff <GFF_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	optional:
	--cutoff <MIN_PROPORTION_OF_TRANSCRIP_COVERED>[90]
	--sample <NAME_OF_EACH_SAMPLE>
	--samtools <SAMTOOLS_PATH>[samtools]
	--bedtools <BED_TOOLS_PATH>[genomeCoverageBed]
	--mincov <MIN_COVERAGE_TO_CONSIDER_POSITION>[1]
	
	bug reports and feature requests: b.pucker@tu-bs.de
					"""

import sys, os, subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- end of imports --- #

def construct_cov_file( bam_file, cov_file, samtools, bedtools, bam_sorted ):
	"""! @brief construct a coverage file """
	
	if bam_sorted:
		sorted_bam_file = bam_file
	else:
		sys.stdout.write( "sorting BAM file ...\n")
		sys.stdout.flush()
		sorted_bam_file = cov_file.split('.cov')[0] + ".sorted.bam"
		cmd = samtools + " sort -m 5000000000 --threads 8 " + bam_file + " > " + sorted_bam_file
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	
	# --- calculate read coverage depth per position --- #
	sys.stdout.write( "calculating coverage per position ....\n" )
	sys.stdout.flush()
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + cov_file
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()


def load_cov_from_file( cov_file ):
	"""! @brief load content of coverage file into dictionary """
	
	coverage = {}	#sequences are keys and coverage values are stored in lists of values
	with open( cov_file, "r" ) as f:
		first_line = f.readline().strip().split('\t')
		prev_seq = first_line[0]
		pos_counter = 1
		covs = []
		while pos_counter < int( first_line[1] ):
			covs.append( 0 )
			pos_counter += 1
		covs.append( float( first_line[1] ) )
		pos_counter += 1
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_seq:	#checks if moving from one sequence to next
				coverage.update( { prev_seq: covs } )
				prev_seq = parts[0] + ""
				covs = []
				pos_counter = 1
				while pos_counter < float( parts[1] ):
					covs.append( 0 )
					pos_counter += 1
				covs.append( float( parts[2] ) )
				pos_counter += 1
			else:
				while pos_counter < float( parts[1] ):
					covs.append( 0 )
					pos_counter += 1
				covs.append( float( parts[2] ) )
				pos_counter += 1
			line = f.readline()
		coverage.update( { prev_seq: covs } )	#adds the coverage values of last sequence
		
	return coverage


def load_transcript_structures_from_gff( gff_file, id_tag="ID" ):
	"""! @brief load exon ranges from given GFF file """
	
	transcript_structures = {}	#dictionary with transcripts as key; dictionaries as values: exon borders + orientation + chromosome name
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript" ]:
					ID = parts[-1].split( id_tag + "=" )[1]
					if ";" in ID:
						ID = ID.split(';')[0]
					transcript_structures.update( { ID: { 'pos': [], 'chr': parts[0], 'orientation': parts[6] } } )
				elif parts[2] == "exon":
					parent = parts[-1].split( "Parent=" )[1]
					if ";" in parent:
						parent = parent.split(';')[0]
					try:
						transcript_structures[ parent ]['pos'].append( [ int( parts[3] ), int( parts[4] ) ] )
					except KeyError:
						sys.stdout.write( "EXON-ERROR: " + line )
						sys.stdout.flush()
						pass
			line = f.readline()
	return transcript_structures


def get_cov_values_per_transcript( coverage, transcript_structures ):
	"""! @brief collect all coverage values per transcript in one list """
	
	covs_per_transcript = {}
	for trans in list( transcript_structures.keys() ):
		try:
			# --- append coverage values of all exon positions to one list --- #
			cov_collection = []
			chromosome_coverage = coverage[ transcript_structures[ trans ][ 'chr' ] ]
			for exon in transcript_structures[ trans ][ 'pos' ]:
				for i in range( exon[0]-1, exon[1] ):	#generates an index list for one exon
					cov_collection.append( chromosome_coverage[ i ] )
		
			# --- check orientation and flip exon order if necessary --- #
			if transcript_structures[ trans ][ 'orientation' ] == '-':
				cov_collection = cov_collection[::-1]
			covs_per_transcript.update( { trans: cov_collection } )
		except:
			pass
	return covs_per_transcript


def write_cov_to_file( cov_input, out_file ):
	"""! @brief write coverage values per transcript into output file """
	
	with open( out_file, "w" ) as out:
		for key in list( sorted( cov_input.keys() ) ):
			values = cov_input[ key ]
			out.write( key + "\t" + ",".join( list( map( str, values ) ) ) + "\n" )


def summarize_across_transcripts( coverages_per_transcript, summary_data_output_file, summary_fig_outout_file, mincov, cutoff ):
	"""! @brief summarize data per transcript """
	
	proportion_of_transcript_covered = {}
	supported_transcripts = []
	for transcript in list( coverages_per_transcript.keys() ):
		counter = 0
		for pos in coverages_per_transcript[ transcript ]:
			if pos >= mincov:
				counter += 1
		proportion_covered = 100.0 * counter / len( coverages_per_transcript[ transcript ] )
		proportion_of_transcript_covered.update( { transcript: proportion_covered } )
		if proportion_covered >= cutoff:
			supported_transcripts.append( transcript )
	
	# --- write into output file --- #
	with open( summary_data_output_file, "w" ) as out:
		out.write( "TranscriptID\tPercentageCovered\n" )
		for transcript in list( sorted( proportion_of_transcript_covered.keys() ) ):
			out.write( transcript + "\t" + str( proportion_of_transcript_covered[ transcript ] ) + "\n" )
	
	# --- generate summary  figure --- #
	values_to_plot = proportion_of_transcript_covered.values()
	fig, ax = plt.subplots()
	ax.hist( values_to_plot )
	ax.set_xlabel( "percentage of transcript covered by RNA-seq" )
	ax.set_ylabel( "number of transcripts" )
	fig.savefig( summary_fig_outout_file, dpi=300 )
	
	return proportion_of_transcript_covered, supported_transcripts


def generate_comparative_plot( collected_data, final_fig_file ):
	"""! @brief generate comparative plot """
	
	fig, ax = plt.subplots()
	
	x_values = list( range( 0, 101 ) )
	for key in list( sorted( collected_data.keys() ) ):
		y_values = []
		for x in x_values:
			y_values.append( 0 )
		for each in collected_data[ key ].values():
			y_values[ int( each ) -1  ] += 1
		ax.plot( x_values, y_values, label=key )
		print( y_values )
	
	ax.legend()
	ax.set_xlabel( "percentage of transcript covered by RNA-seq" )
	ax.set_ylabel( "number of transcripts" )
	
	fig.savefig( final_fig_file, dpi=300 )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	
	if not output_folder[-1] == "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--samtools' in arguments:
		samtools = arguments[ arguments.index('--samtools')+1 ]
	else:
		samtools = "samtools"
	
	if '--bedtools' in arguments:
		bedtools = arguments[ arguments.index('--bedtools')+1 ]
	else:
		bedtools = "genomeCoverageBed"
	
	if '--bam_is_sorted' in arguments:
		bam_sorted = True
	else:
		bam_sorted = False
	
	if '--cutoff' in arguments:
		cutoff = int( arguments[ arguments.index('--cutoff')+1 ] )
	else:
		cutoff = 90
	
	if '--mincov' in arguments: #how many bases of a transcript were sequenced in total (sum of coverage of individual positions)
		mincov = int( arguments[ arguments.index('--mincov')+1 ] )
	else:
		mincov=1
	
	if '--sample' in arguments:
		sample = arguments[ arguments.index('--sample')+1 ]
		if "," in sample:
			samples = sample.split(',')
	else:
		if '--bam' in arguments:
			number = arguments[ arguments.index('--bam')+1 ].count(',')+1
			samples = list( map( str, range( number ) ) )
		else:
			number = arguments[ arguments.index('--cov')+1 ].count(',')+1
			samples = list( map( str, range( number ) ) )
	
	if '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]	#path to BAM input file
		#generate coverage file / also allow start from coverage file
		if "," in bam_file:
			bam_files = bam_file.split(',')
		else:
			bam_files = [ bam_file ]
		cov_files = []
		for idx, bam_file in enumerate( bam_files ):
			cov_file = output_folder + samples[idx] + ".cov"
			if not os.path.isfile( cov_file ):
				sys.stdout.write( "constructing coverage file "+ samples[ idx ] +" ...\n" )
				sys.stdout.flush()
				construct_cov_file( bam_file, cov_file, samtools, bedtools, bam_sorted  )
				sys.stdout.write( "...done\n" )
				sys.stdout.flush()
			cov_files.append( cov_file )
	else:
		cov_file = arguments[ arguments.index('--cov')+1 ]
		if "," in cov_file:
			cov_files = cov_file.split(',')
		else:
			cov_files = [ cov_file ]
	

	# --- load gene structures from GFF (exon ranges of representative transcript) --- #
	sys.stdout.write( "loading transcript structures ...\n" )
	sys.stdout.flush()
	transcript_structures = load_transcript_structures_from_gff( gff_file )
	sys.stdout.write( "number of identified transcripts: "+ str( len( list( transcript_structures.keys() ) ) ) +"\n" )
	sys.stdout.write( "...done\n" )
	sys.stdout.flush()
	
	# --- iteration over all samples --- #
	collected_data = {}
	for oidx, cov_file in enumerate( cov_files ):	#iterate over all coverage files that have been created in the previous step
		# --- load coverage from file --- #
		sys.stdout.write( "loading coverage from COV file "+ samples[ oidx ] +" ...\n" )
		sys.stdout.flush()
		coverage = load_cov_from_file( cov_file )
		sys.stdout.write( "...done\n" )
		sys.stdout.flush()
		
		#collect coverage per position per transcript
		sys.stdout.write( "collecting coverage per transcript "+ samples[ oidx ] +" ...\n" )
		sys.stdout.flush()
		coverages_per_transcript = get_cov_values_per_transcript( coverage, transcript_structures )
		sys.stdout.write( "...done\n" )
		sys.stdout.flush()
		
		# --- write coverage per transcript into output file --- #
		cov_per_transcript_out_file = output_folder + samples[ oidx ] + ".cov_per_transcript.txt"
		write_cov_to_file( coverages_per_transcript, cov_per_transcript_out_file )
		
		#summarize the covered proportions
		summary_data_output_file = output_folder + samples[ oidx ] + ".summary.txt"
		summary_fig_outout_file = output_folder + samples[ oidx ] + ".summary.png"
		summary, supported_transcripts = summarize_across_transcripts( coverages_per_transcript, summary_data_output_file, summary_fig_outout_file, mincov, cutoff )
		
		#supported_transcripts = list of IDs to keep
		
		
		# --- store summary information --- #
		collected_data.update( { samples[ oidx ]: summary } )
	
	# --- generate final comparative figure --- #
	final_fig_file = output_folder + "comparative_plot.png"
	generate_comparative_plot( collected_data, final_fig_file )


if '--bam' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
