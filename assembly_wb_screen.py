### Boas Pucker ###
### b.pucker@tu-bs.de ###

### some functions taken from KIPEs and MYB_annotator ###

__version__ = "v0.1"

__reference__ = "Pucker, 2023"

__usage__ = """
					Screen genome assembly for contamination against white and black list """ + __version__ + """("""+ __reference__ +""")
					
					Usage:
					python3 assembly_wb_screen.py
					--in <ASSEMBLY_FILE>
					--white <WHITE_LIST_FASTA_FILE>
					--black <BLAST_LIST_FASTA_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--tmp <TMP_FOLDER>[output folder]
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, sys, re, glob, subprocess

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load sequences of given FASTA file into dictionary with sequence IDs as keys and sequences as values """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:	#take only the space-free part of the header (if space present)
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def generate_fragment_file( assembly_file, fragment_file, fragment_size ):
	"""! @brief generate file with sequence fragments based on assembly file """
	
	assembly = load_sequences( assembly_file )
	
	with open( fragment_file, "w" ) as out:
		for key in list( assembly.keys() ):
			chunks = [ assembly[ key ][ i:i+fragment_size ] for i in range( 0, len( assembly[ key ] ), fragment_size ) ]
			for idx, chunk in enumerate( chunks ):
				out.write( '>' + key + "_%_" + ( str( idx ).zfill(4) ) + "\n" + chunk + "\n" )


def load_BLAST_results( blast_result_file ):
	"""! @brief load best BLAST hit per query from given BLAST result file"""
	
	hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > hits[ parts[0] ]['score']:
					hits[ parts[0] ] = { 'score': float( parts[-1] ), 'sim': float( parts[2] ), 'len': int( parts[3] ) }
			except KeyError:
				hits.update( { parts[0]: { 'score': float( parts[-1] ), 'sim': float( parts[2] ), 'len': int( parts[3] ) } } )
			line = f.readline()
	return hits


def generate_BLAST_result_output_file( white_hits, black_hits, headers, detail_output_file, cutoff_ratio ):
	"""! @brief generate a summary table with BLAST hits against black and white """
	
	with open( detail_output_file, "w" ) as out:
		out.write( "\t".join( [ "ContigPartID", "WhiteScore", "WhiteSim", "WhiteLen", "BlackScore", "BlackSim", "BlackLen", "ScoreRatio", "Status" ] ) + "\n" )
		for header in headers:
			new_line = [ header ]
			
			# --- white --- #
			try:
				white = white_hits[ header ]
				new_line.append( str( white['score'] ) )
				new_line.append( str( white['sim'] ) )
				new_line.append( str( white['len'] ) )
				white_score = white['score']
			except KeyError:
				new_line.append( "0" )
				new_line.append( "0" )
				new_line.append( "0" )
				white_score = 1
			
			# --- black --- #
			try:
				black = black_hits[ header ]
				new_line.append( str( black['score'] ) )
				new_line.append( str( black['sim'] ) )
				new_line.append( str( black['len'] ) )
				black_score = black['score']
			except KeyError:
				new_line.append( "0" )
				new_line.append( "0" )
				new_line.append( "0" )
				black_score = 1
			
			# --- compare best white against best black hit to perform contamination classification --- #
			new_line.append( str( white_score/black_score ) )
			if white_score/black_score > cutoff_ratio:
				new_line.append( "OK" )
			elif white_score/black_score >= 1:
				new_line.append( "UNCLEAR" )
			else:
				new_line.append( "WARNING" )
			
			out.write( "\t".join( new_line )+"\n" )


def main( arguments ):
	"""! @brief run everything """
	
	assembly_file = arguments[ arguments.index('--in')+1 ]
	white_file = arguments[ arguments.index('--white')+1 ]
	black_file = arguments[ arguments.index('--black')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	
	if '--tmp' in arguments:
		tmp_folder = arguments[ arguments.index('--tmp')+1 ]
	else:
		tmp_folder = output_folder + "tmp/"
	
	if tmp_folder[-1] != "/":
		tmp_folder += "/"
	
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )
	
	
	fragment_size = 10000	#10kb (size of individual fragments for BLAST search)
	evalue = 0.0001
	cpus = 1
	cutoff_ratio = 2	#score ratio of best white vs. best black hit to consider sequence as clean
	
	
	fragment_file = output_folder + "fragments.fasta"
	if not os.path.isfile( fragment_file ):
		generate_fragment_file( assembly_file, fragment_file, fragment_size )	#split assembly into parts
	
	# --- run BLAST vs. white --- #
	white_blast_result_file = tmp_folder + "white_blast_result_file.txt"
	if not os.path.isfile( white_blast_result_file ):
		white_blast_db = tmp_folder + "white_db"
		p = subprocess.Popen( args= "makeblastdb -in " + white_file + " -out " + white_blast_db + " -dbtype nucl", shell=True )
		p.communicate()
		
		p = subprocess.Popen( args= "blastn -query " + fragment_file + " -db " + white_blast_db + " -out " + white_blast_result_file + " -outfmt 6 -evalue " + str( evalue ) + " -num_threads " + str( cpus ), shell=True )
		p.communicate()
	
	# --- run BLAST vs. black --- #
	black_blast_result_file = tmp_folder + "black_blast_result_file.txt"
	if not os.path.isfile( black_blast_result_file ):
		black_blast_db = tmp_folder + "black_db"
		p = subprocess.Popen( args= "makeblastdb -in " + black_file + " -out " + black_blast_db + " -dbtype nucl", shell=True )
		p.communicate()
		
		p = subprocess.Popen( args= "blastn -query " + fragment_file + " -db " + black_blast_db + " -out " + black_blast_result_file + " -outfmt 6 -evalue " + str( evalue ) + " -num_threads " + str( cpus ), shell=True )
		p.communicate()

	# --- load BLAST results --- #
	white_hits = load_BLAST_results( white_blast_result_file )
	black_hits = load_BLAST_results( black_blast_result_file )
	
	# --- analyze results --- #
	headers = sorted( list( load_sequences( fragment_file ).keys() ) )
	detail_output_file = output_folder + "BLAST_result_details.txt"
	generate_BLAST_result_output_file( white_hits, black_hits, headers, detail_output_file, cutoff_ratio )



if '--in' in sys.argv and '--out' in sys.argv and '--white' in sys.argv and '--black' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
