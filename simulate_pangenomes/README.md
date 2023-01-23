# Fragmentation and Incompleteness simulations

From manuscript: 
```
Critical assessment of pan-genomic analysis of metagenome-assembled genomes 
Tang Li, Yanbin Yin
Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac413, https://doi.org/10.1093/bib/bbac413
Published: 17 September 2022
```

## Process:

Complete genomes need to be fragmented and then artificially have parts of the broken contigs removed: 

```
for f in ../fasta/*; do python fragmentation.py $f; done
for f in ../fasta/*cut; do python incompleteness.py $f ; done 
```

Fragmentation and incompleteness are set in the scripts (currently set at 100 fragments, 80% completion). 
