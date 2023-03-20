# Benga
Bacterial Epidemiology NGs Analysis (BENGA) framework and pipeline.

# Requirements
* python 3.6+
  * biopython
  * scipy
  * numba
  * numpy
  * pandas
  * fastcluster
  * matplotlib
* ncbi-blast
* prodigal
* prokka
* roary

# Usage
### ***cgMLST profiling***
* You can create yourself training file by prodigal. 
```
profiling -i fasta_file -o result.tsv -d <database> --prodigaltf <prodigal training file>
```
### ***Pan genome analysis***
```
create_scheme.py -i input_path -o output_path
```
### ***Defining the cgMLST scheme***
* Defind cgMLST scheme by gene occurrence, default is 95%.
```
extract_scheme.py -i <output of create_scheme.py> -o <path of database>
```
### ***Draw dendrogram***
```
cluster.py -i <profile1> <profile2> <profile3> ... -o <output_path>
```
