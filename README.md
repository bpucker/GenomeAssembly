# GenomeAssembly
collection of tools for the analysis and cleaning of genome assemblies




## FASTQ statistics calculation

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

`--in` specifies a BAM file containing folder.


## Assembly statistics calculation

```
Usage:
  python3 contig_stats3.py --in <DIR> --out <DIR>
  
  --in   STR   BAM input folder
  --out  STR   BAM output file
```

`--in` specifies a BAM file containing folder.




## Screen assembly for contamination (white list and black list)

```
Usage:
  python3 assembly_wb_screen.py --in <DIR> --out <DIR>
  
  --in   STR   BAM input folder
  --out  STR   BAM output file
```

`--in` specifies a BAM file containing folder.



# References

This repository.

