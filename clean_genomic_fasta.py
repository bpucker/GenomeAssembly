### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

__version__ = "v0.25"

__usage__ = """
					python3 clean_genomic_fasta.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == ">":
					if " " in line:
						tmp = line.split(' ')[0].replace("|", "_")
						if chr(9) in tmp:
							tmp = tmp.split( chr(9) )[0]
						out.write(  tmp + "\n" )
					else:
						if chr(9) in line:
							out.write( line.split( chr(9) )[0].replace("|", "_") + "\n" )
						else:
							out.write( line.strip().replace("|", "_") + "\n" )
				else:
					out.write( line.strip() + "\n" )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
