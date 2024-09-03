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
    #filter for highly variable genes