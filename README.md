# Insertion and deletion Reconstruction tool


This function reconstructs insertions on each branch in RSV phylogenetic trees, 
compares them to deletions and outputs a CSV file comparing insertions and deletions.


### Input Data

The function inputs include:
 
 * genbank format genome reference
 
 * annotated tree (json format)
 
 * tree (nwk format)
 
 * insertion file (csv format)
 
 * aligned G amino acid file (fasta)
 
 The reference files are available in the data/ folder. The rest of the inputs can be generated by running the [without\_G_workflow](https://github.com/LauraU123/without_G_workflow)
 

### Output Data

The output consists of a TSV file with columns for insertion and deletion length and location in the G gene. 


### Running the function

This function can be run manually the command line as:

```
	python3 insertions_and_deletions.py \

	--tree {tree file nwk format} \

	--json {annotated tree json format} \

	--reference {annotated genome reference gbk format} \

	--insertions {insertions file CSV format} \
	
	--alignedgene {aligned G amino acid file} \

	--output {output file TSV format}
```

It can also be run using Snakemake with
```
	snakemake --cores all
```
