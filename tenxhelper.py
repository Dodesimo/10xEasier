from flask import request
import os
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt

def download():
    download = request.files.getlist("uploadbutton")

    #save file, return file name
    fileName = os.path.join("/Users/devammondal/Downloads/Appseq/App/downloads", download.filename)
    download.save(fileName)
    return fileName

def generateAnn(fileName):
    #create data object, and return it.
    data = sc.read_10x_mtx('downloads', var_names='gene_symbols', cache=True)
    return data

def preprocess(data):
    #minimum gene and cell processing for data.
    sc.pp.filter_cells(data, min_genes=200)
    sc.pp.filter_genes(data, min_cells=3)

    #process mitochondrial
    data.var['mt'] = data.var_names.str.startswith('MT-')

def qualityControl(data):
    #introduce new column for quality control
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


def generateGraphics(data):
    #generate three figures (one violin, second assessing mitochondrial counts, third assessing genes per count), and save those in the outputs folder.
    figure = sc.pl.violin(data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    figure.savefig(os.path.join("/Users/devammondal/Downloads/Appseq/App/outputs", "qualityControlViolin"))
    figure = sc.pl.scatter(data, x='total_counts', y='pct_counts_mt')
    figure.savefig(os.path.join("/Users/devammondal/Downloads/Appseq/App/outputs", "qualityControlMitochondrial"))
    figure = sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts')
    figure.savefig(os.path.join("/Users/devammondal/Downloads/Appseq/App/outputs", "qualityControlGeneCount"))

def normalization(data):
    #Do log count normalization
    sc.pp.log1p(data)

def clustering(data):
    #establish criteria for highly variable genes
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)

    #filter out the highly variable genes
    data = data[:, data.var.highly_variable]

    #regressing out total_counts, pct_counts_mt (removing correlating effects between them)
    sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])

    #scale results
    sc.pp.scale(data, max_value=10)

    #PCA analysis, single value decomposition using arpack.
    sc.tl.pca(data, svd_solver='arpack')

    #Calculate variance ratio figure of PCAs.
    figure = sc.pl.pca_variance_ratio(data, log=True)
    figure.savefig(os.path.join("/Users/devammondal/Downloads/Appseq/App/outputs", "principalComponentsRatios"))

    #Cluster.
    sc.pp.neighbors(data, n_neighbors=10, n_pcs=20)

    #Create both leiden and zika colorings of umap configuration
    figure = sc.pl.umap(data, color='leiden')
    figure.savefig(os.path.join("/Users/devammondal/Downloads/Appseq/App/outputs", "uMapLeidenClustering"))



