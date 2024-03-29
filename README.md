[![DOI](https://zenodo.org/badge/623857647.svg)](https://zenodo.org/badge/latestdoi/623857647)

# GenomeAssembly
collection of tools for the analysis and cleaning of genome assemblies




## FASTQ statistics calculation
This script calculates statistics of a FASTQ file including the number of reads, the average read length, the total number of sequenced bases and the read length N50. Additionally, a read length histogram and a quality vs. read length plot can be generated.

```
Usage:
  python3 FASTQ_stats3.py --in <FILE> | --in_dir <DIR>
  
  mandatory:
  --in      STR   Input FASTQ file
  --in_dir  STR   Input folder
  
  optional:
  --rfig    STR   Read length histogram figure filename
  --cutoff  STR   Read length cutoff for histogram (kb) [100]
	
  --qfig    STR   Quality vs. read length figure filename
  --lencut  STR   Upper read length cutoff (kb) [200]
  --qualcut STR   Upper quality cutoff (phred) [40]
```

`--in` specifies a FASTQ input file that will be analyzed. The file should be gzip compressed.

`--in_dir` specifies a FASTQ file containing input folder. Each (gzip compressed) FASTQ file in the folder will be analyzed. Supported file extensions: .fq, .fastq, .fq.gzip, fq.gz, fastq.gzip, .FQ, .FASTQ, .FQ.GZIP, .FASTQ.GZIP, .FQ.GZ, and .FASTQ.GZ.

`--rfig` specifies the filename of a read length histogram figure. Inclusion of this argument triggers the generation of this figure. Defaul: off.

`--cutoff` specifies the upper read length cutoff of the read length histogram figure. Default: 100 (kb).

`--qfig` specifies filename of a quality vs. read length figure. Inclusion of this argument triggers the generation of this figure. Defaul: off.

`--lencut` specifies an upper read length cutoff for the qality vs. read length figure. Default: 200 (kb).

`--qualcut` specifies an upper quality cutoff for the quality vs. read length figure. Default: 40.



## Assembly statistics calculation
This script calculates the statistics of a given assembly FASTA file. This includes trimming short contigs and cleaning the contig names. Statistics include the number of contigs, the N50, N90, and total assembly size. If an expression file is provided, the calculation of Ex90N50 can be performed.

```
Usage:
  python3 contig_stats3.py --in <FILE>
  
  --in   STR   Input FASTA file
  
  optional:
  --min  INT   Minimal contig length [1000]
  --out  STR   Output folder
  --exp  STR   Expression file
```

`--in` specifies a FASTA file that will be analyzed. A trimmed FASTA file and a statistics file will be placed next to the input file.

`--min` specifies the minimal contig length. Default: 1000 (bp).

`--out` specifies an output folder. The trimmed FASTA file and the statistics file will be placed in this folder. Default: off.

`--exp` specifies an expression file (normalized expression). Default: none.


## Clean genomic FASTA
This script cleans the header names of a given FASTA file by splitting at the first white space chracter (space or tab). Some special characters are also replaced by underscores.

```
Usage:
  python3 clean_genomic_fasta.py --in <FILE> --out <FILE>
  
  --in   STR   Input FASTA file
  --out  STR   Output FASTA file
```

`--in` specifies a FASTA file that will be processed.

`--out` specifies the output FASTA file.




## Screen assembly for contamination (white list and black list)
This script screens an assembly for contaminations. All contigs are separted into small blocks and compared against a white list and black list of sequences. This allows the identification of contamination contigs or chimeric contigs that contain contaminations.


```
Usage:
  python3 assembly_wb_screen.py --in <FILE> --out <DIR> --white <FILE> --black <FILE>
  
  --in     STR   Input FASTA folder
  --out    STR   BAM output file
  --white  STR   White list FASTA file
  --black  STR   Black list FASTA file
  
  optional:
  --tmp    STR   Temp output folder
```

`--in` specifies an assembly FASTA file that will be screened for contaminations.

`--out` specifies an output folder where all result files can be placed. The folder will be created if it does not exist already.

`--white` specifies a FASTA file with white sequences that will be used to identify valid contigs and contig parts.

`--black` specifies a FASTA file with black sequences that will be used to identify invalid (contamination) contigs and contig parts.

`--tmp` specifies a temporary output folder. The `--out` folder is used for temporary files if this argument is not used.



## Identify best supported gene models based on RNA-seq coverage ##
This script analyzes the RNA-seq coverage of predicted transcripts. The results allow conclusions about the valididty of predicted transcripts.


```
Usage1
python3 RNAseq_cov_analysis.py --bam <FILE> --gff <FILE> --out <DIR>

Usage2
python3 RNAseq_cov_analysis.py --cov <FILE> --gff <FILE> --out <DIR>

Mandatory:
--bam        STR   BAM input file    |  --cov   STR    Coverage input file
--gff        STR   GFF input file
--out        STR   Output folder

Optional:
--cutoff     INT    Percentage of transcript covered[90]
--samtools   STR    Path to samtools [samtools]
--bedtools   STR    Path to bedtools [bedtools]
--mincov     INT    Minimal coverage per position [1]
```

`--bam` specifies a BAM input file. This will be converted into a coverage file (COV) to analyze the distribution of reads across transcripts. This argument can also be used to provide a comma-separated list of files for automatic processing of large batches of files.

`--cov` specifies a coverage file (COV) to analyze the distribution of reads across transcripts. This argument can also be used to provide a comma-separated list of files for automatic processing of large batches of files.

`--gff` specifies a GFF file that contains the exon positions for the RNA-seq coverage analysis.

`--out` specifies an output folder. If this folder does not exist, it will be created.

`--cutoff` specifies a minimal proportion of the transcript that needs to be covered by RNA-seq. Default value is 90.

`--samtools` specifies the full path to samtools. Default: samtools.

`--bedtools` specifies the full path to bedtools. Default: bedtools.

`--mincov` specifies the minimal coverage of a position in order to consider it. Default: 1.



# References

Meckoni S. N., Nass B., Pucker B. (2023). Phylogenetic placement of _Ceratophyllum submersum_ based on a complete plastome sequence derived from nanopore long read sequencing data. bioRxiv 2023.06.27.546741; doi: [10.1101/2023.06.27.546741](https://doi.org/10.1101/2023.06.27.546741).


