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
  * click
* [ncbi-blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
* [prodigal](https://github.com/hyattpd/Prodigal)
* [prokka](https://github.com/tseemann/prokka)
* [roary](https://github.com/sanger-pathogens/Roary)

# Usage
### ***cgMLST profiling***
Convert genome sequence to cgMLST profile. We provide some organism database in folder **scheme**. 
```
Benga.py cgmlst_profiling -i <genome.fasta> -o <profile.tsv> -d <organism.db>
```
It is possible to use the same gene prediction model for each genome. This will improve the consistency of gene prediction but any biases present in the model will be present in all genome annotations.  
If you want to do this, you can add `--prodigaltf <prodigal training file>`, we also provide some model in the folder **prodigaltf**, or you can make it with `prodigal`
### ***Create wgMLST scheme***
If you want to create your cgMLST scheme, you have to create a wgMLST scheme in first. Given a directory containing genome assemblies in FASTA format to create it.
```
Benga.py create_scheme -i <path/of/input> -o <path/of/output>
```
### ***Defining the cgMLST scheme***
Define cgMLST scheme by gene occurrence, default is 95%. It will make a sqlite database for `cgmlst_profiling`.
```
Benga.py extract_scheme -i <output/of/create_scheme.py> -o <prefix/of/database>
```
### ***Draw dendrogram***
Fast cluster genome that use single linkage. Fast cluster that use single linkage. The cgMLST profiles must in input directory.
```
Benga.py cgmlst_cluster -i <path/of/input> ... -o <path/of/output>
```

# Citation
If you use Benga please cite:  

**Chen YS, Tu YH, Chen BH, Liu YY, Hong YP, Teng RH, Wang YW, Chiou CS. cgMLST@Taiwan: A web service platform for Vibrio cholerae cgMLST profiling and global strain tracking. J Microbiol Immunol Infect. 2022 Feb;55(1):102-106. doi: 10.1016/j.jmii.2020.12.007. Epub 2021 Jan 15. PMID: 33485793.**
