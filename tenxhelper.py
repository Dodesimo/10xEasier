from flask import request
import os
import pandas as pd
import scanpy as sc
import matplotlib
from matplotlib import pyplot as plt

matplotlib.use('Agg')


def download():
    files = request.files.getlist("uploadbutton")
    saved_files = []

    for file in files:
        fileName = os.path.join("/Users/devammondal/Downloads/Appseq/App/downloads", file.filename)
        file.save(fileName)
        saved_files.append(fileName)

    return saved_files


def generateAnn():
    # create data object, and return it.
    data = sc.read_10x_mtx('downloads', var_names='gene_symbols', cache=True)
    return data


def preprocess(data):
    # minimum gene and cell processing for data.
    sc.pp.filter_cells(data, min_genes=200)
    sc.pp.filter_genes(data, min_cells=3)

    # process mitochondrial
    data.var['mt'] = data.var_names.str.startswith('MT-')
    return data


def qualityControl(data):
    # introduce new column for quality control
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return data


def generateGraphics(data):
    output_dir = "/Users/devammondal/Downloads/Appseq/App/downloads"
    # Violin plot for quality control metrics
    sc.pl.violin(data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    plt.savefig(os.path.join(output_dir, "qualityControlViolin.png"))
    plt.close()  # Close the figure to avoid overlap with other plots

    # Scatter plot for mitochondrial percentage vs total counts
    sc.pl.scatter(data, x='total_counts', y='pct_counts_mt')
    plt.savefig(os.path.join(output_dir, "qualityControlMitochondrial.png"))
    plt.close()

    # Scatter plot for gene counts vs total counts
    sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts')
    plt.savefig(os.path.join(output_dir, "qualityControlGeneCount.png"))
    plt.close()


def normalization(data):
    # Do log count normalization
    sc.pp.log1p(data)
    return data


def clustering(data):
    output_dir = "/Users/devammondal/Downloads/Appseq/App/downloads"

    # establish criteria for highly variable genes
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # filter out the highly variable genes
    data = data[:, data.var.highly_variable]

    # regressing out total_counts, pct_counts_mt (removing correlating effects between them)
    sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])

    # scale results
    sc.pp.scale(data, max_value=10)

    # PCA analysis, single value decomposition using arpack.
    sc.tl.pca(data, svd_solver='arpack')

    # Calculate variance ratio figure of PCAs.
    sc.pl.pca_variance_ratio(data, log=True)
    plt.savefig(os.path.join(output_dir, "pcavarratio.png"))
    plt.close()

    # Cluster.
    sc.pp.neighbors(data, n_neighbors=10, n_pcs=20)
    sc.tl.umap(data)
    sc.pl.umap(data)
    plt.savefig(os.path.join(output_dir, "umap.png"))
    plt.close()

    # Create both leiden and zika colorings of umap configuration
    sc.tl.leiden(data, resolution=0.25)
    sc.pl.umap(data, color=['leiden'])
    plt.savefig(os.path.join(output_dir, "umapleiden.png"))
    plt.close()


