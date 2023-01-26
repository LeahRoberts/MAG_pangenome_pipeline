"""
script to make gene curve plots from pangenome output
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys

pangenome_data = sys.argv[1]

data = pd.read_csv(pangenome_data, sep='\t', header=None)
dataDf = pd.DataFrame(data)
gene_proportionsDf = dataDf.iloc[:, [2]]
plt.hist(gene_proportionsDf, bins=10)
plt.savefig("pangenome_curve.png")
#plt.show()

