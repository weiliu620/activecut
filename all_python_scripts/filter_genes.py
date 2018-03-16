import sys
import os
import numpy as np
import csv

def make_gene_filter(pacall_file):
    """
    filter genes by expression.

    Some expressions are overly expressed, some others are not expressed at all. filter each row of the data matrix (genes x samples) so only moderately expressed remains.
    input: the PACall.csv, downloaed from the Allen institute website.
    output: a boolean np.ndarray that have True for included genes, and False for exclued.
    
    """
    
    pacall = np.genfromtxt(pacall_file, delimiter = ',', dtype = 'int')
    pacall = pacall[:,1:]
    rsum = pacall.sum(axis = 1)
    filter = (rsum < 940) & (rsum > 10)
    return np.array(filter)

def make_gray_filter(ontology_file, sample_anno_file, mask_id):
    """
    filter sample columns by anatomical structure ('4006' for gray matter).
    """
    
    # create a dictionary of if a region is gray matter.
    gray_dic = {}
    with open(ontology_file, 'rb') as f:
        reader = csv.reader(f)
        reader.next() # skip header
        for row in reader:
            if mask_id in row[6]:
                gray_dic[row[0]] = True
            else:
                gray_dic[row[0]] = False
    
    # filter the regions.
    filter = list()
    sample_coord = list()
    with open(sample_anno_file, 'rb') as f:
        reader = csv.reader(f)
        reader.next()
        for row in reader:
            filter.append(gray_dic[row[0]])
            sample_coord.append([row[10], row[11], row[12]])

    
    return np.array(filter)

def filter_genes(filter, expression_file, out_file):
    """
    apply the filter to gene expression matrix.
    """
    
    expressions = np.genfromtxt(expression_file, delimiter = ',')
    expressions_new = expressions[filter, :]
    np.savetxt(out_file, expressions_new, delimiter = ',', fmt = "%.6f")

def filter_coord(sample_anno_file, filter):
    """
    filtering the sample annotation file by the anatomical filter.

    Output a list of coordinates.
    """
    all_coord = np.genfromtxt(sample_anno_file, skip_header = 1, delimiter = ',', usecols = (10, 11, 12))
    # apply the filter
    sample_coord_new = all_coord[filter,:]
    
    return sample_coord_new

def filter_samples(filter, expression_file, out_file):
    """
    Apply filter to the sample column of gene expression array.
    """
    expressions = np.genfromtxt(expression_file, delimiter = ',')

    filter = np.concatenate((np.array([True]), filter))
    expressions_new = expressions[:,filter]
    np.savetxt(out_file, expressions_new, delimiter = ',', fmt = "%.6f")

