# Insertion and deletion Reconstruction tool


This function reconstructs insertions on each branch in RSV phylogenetic trees, 
compares them to deletions and outputs a CSV file comparing insertions and deletions.

### Input Data

The function inputs include:
 
 * genbank format genome reference
 
 * annotated tree (json format)
 
 * tree (nwk format)
 
 * insertion file (csv format)

### Output Data

The output consists of a TSV file with columns for insertion and deletion length and location in the G gene. 


### Running the function

This function can be run from the command line as:

```
	python3 insertions_and_deletions.py \

	--tree {tree file nwk format} \

	--json {annotated tree json format} \

	--reference {annotated genome reference gbk format} \

	--insertions {insertions file CSV format} \

	--output {output file TSV format}
```
