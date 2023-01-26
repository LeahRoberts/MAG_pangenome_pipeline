"""
script to summarise pangenomes
"""

import sys
import pandas as pd

core_perc = float(sys.argv[1])
gene_desc = sys.argv[2]
matrix = sys.argv[3]
outputfile = sys.argv[4]

# import matrix as pandas df
gene_matrixDf = pd.read_csv(matrix, sep=',', header=0, index_col=0)
genes = (list(gene_matrixDf)[0:])
MAG_list = gene_matrixDf.index.to_list()

# put gene descriptions in dictionary:
gene_description_dictionary = {}
with open(gene_desc, "r") as desc_in:
    for line in desc_in:
        gene_name = line.split()[0]
        gene_name = gene_name.replace(">", "")
        description = line.split()[1:]
        description = "_".join(description)
        gene_description_dictionary[gene_name] = description

for gene in genes:
    # sum number of MAGs the gene appears in:
    count = gene_matrixDf[gene].sum()
    prevalence = float(count) / float(len(MAG_list))
    if prevalence >= core_perc:
        pan_type = "core"
    else:
        pan_type = "accessory"
    with open(outputfile, "a") as fout:
        fout.write("%s\t%s\t%s\t%s\n" % (gene, gene_description_dictionary[gene], prevalence, pan_type))


