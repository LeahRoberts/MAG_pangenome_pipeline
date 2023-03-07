# MAG_pangenome_pipeline

A pipeline written in Snakemake to automatically generate pangenomes from metagenome assembled genomes (MAGs). 

## Dependencies: 

* Snakemake (to get Snakemake to work with LSF click [here](https://github.com/Snakemake-Profiles/lsf#install))
* mmseqs2
* Bakta (if running annotation)
* Biopython
* BCBio
* Pandas

**NOTE:** Conda is used to call different environments and dependencies (see Snakemake file).

## Input: 

* Unannotated MAGs (will be annotated with Bakta); OR
* GFF (with fasta)


## Quick start: 

Update `config.yaml` to specify workflow and directory paths. 

Run snakemake (bash script for running on cluster using LSF):

```
bash submit_lsf.sh
```


