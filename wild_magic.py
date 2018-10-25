#!/usr/bin/env python3 
'''
load modules python/3.6.0 and tcltk/8.6.4

To run:
python3 wild_magic.py ~/path/to/inputfile.csv dataset/

'''

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import magic
import sys


if len(sys.argv) != 2:
        print("Usage: python3 wild_magic.py ~/path/to/inputfile.csv dataset/")
        quit()

#if '/' in sys.argv[2]:
#        dataset = sys.argv[2]
#elif '/' not in sys.argv[2]:
#        dataset = sys.argv[2] + '/'

sys.stdout = open('magic.log', 'w')

outfile = os.path.basename(sys.argv[1])
outfile = outfile[:-4]
outfile = outfile + '_'

#outfile = dataset + outfile

print("Running from_csv()")
scdata = magic.mg.SCData.from_csv(os.path.expanduser(sys.argv[1]), data_type='sc-seq', normalize=False)
print("Read from_csv() success")
print("Shape of scdata._data: %d x %d" % (scdata._data.shape[0], scdata._data.shape[1]))


# Set filtering thresholds
#fig, ax = scdata.plot_molecules_per_cell_and_gene()
#data already filtered

# Minimum molecules/cell value
CELL_MIN = 0

# Maximum molecules/cell values
CELL_MAX = 1000000

# Minimum number of nonzero cells/gene 
# (None if no filtering desired)
GENE_NONZERO = None

# Minimum number of molecules/gene
# (None if no filtering desired)
GENE_MOLECULES = None

print("Running filter_scseq_data()")
scdata.filter_scseq_data(filter_cell_min=CELL_MIN,
		         filter_cell_max=CELL_MAX,
                         filter_gene_nonzero=GENE_NONZERO,
		         filter_gene_mols=GENE_MOLECULES)
print("filter_sceq_data() success")
print("Shape of scdata._data: %d x %d" % (scdata._data.shape[0], scdata._data.shape[1]))

# Normalize
print("Running normalize_scseq_data()")
scdata = scdata.normalize_scseq_data()
print("normalize_sc_seq_data() success")
print("Shape of scdata._data: %d x %d" % (scdata._data.shape[0], scdata._data.shape[1]))

print("Running plot_pca_variance_expplained()")
fig, ax = scdata.plot_pca_variance_explained(n_components=40, random=True)
fig.savefig(outfile + 'pca.pdf')
print("plot_pca_variance_expplained() success")

print("Running to_csv() on scdata._data")
dat = scdata._data
dat.to_csv(outfile + "magicdata_before.csv")
print("to_csv() success")

# MAGIC
print("Running run_magic()")
scdata.run_magic(n_pca_components=20, random_pca=True, t=6, k=30, ka=10,epsilon=1, rescale_percent=99)
print("run_magic() success")
print("Shape of scdata._data: %d x %d" % (scdata._data.shape[0], scdata._data.shape[1]))

# save to csv after magic
print("Running to_csv() on scdata._magic._data")
dat = scdata._magic._data
dat.to_csv(outfile + "normalized_magicdat.csv")
print("to_csv() success")

def scatter_genes(gene_one, gene_two, gene_color, magicplot):
    if magicplot == False:
        fig, ax = scdata.scatter_gene_expression([gene_one, gene_two], color=gene_color)
        ax.set_xlabel(gene_one)
        ax.set_ylabel(gene_two)
        fig.savefig(outfile + gene_color + "_before_magic.pdf")
    elif magicplot == True:
        fig, ax = scdata.magic.scatter_gene_expression(["MAGIC " + gene_one, "MAGIC " + gene_two], color="MAGIC " + gene_color)
        ax.set_xlabel("MAGIC " + gene_one)
        ax.set_xlabel("MAGIC " + gene_two)
        fig.savefig(outfile + gene_color + "_after_magic.pdf")
    return

print("Running scatter_genes()")
scatter_genes("VIM", "CDH1", "ZEB1", False)
scatter_genes("VIM", "CDH1", "ZEB1", True)

scatter_genes("VIM", "CDH1", "ZEB2", True)
scatter_genes("VIM", "CDH1", "ZEB2", True)

scatter_genes("VIM", "CDH1", "TWIST1", True)
scatter_genes("VIM", "CDH1", "TWIST1", True)

scatter_genes("VIM", "CDH1", "TWIST2", False)
scatter_genes("VIM", "CDH1", "TWIST2", True)

scatter_genes("VIM", "CDH1", "SNAI1", True)
scatter_genes("VIM", "CDH1", "SNAI1", True)

scatter_genes("VIM", "CDH1", "SNAI2", True)
scatter_genes("VIM", "CDH1", "SNAI2", True)
print("scatter_genes() done")

sys.exit()

