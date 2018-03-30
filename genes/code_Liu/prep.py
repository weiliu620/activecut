import sys
import os
import numpy as np
import csv
import socket
import shutil
import subprocess
if socket.gethostname() == 'wukong':
    import matplotlib.pyplot as plt
    import nibabel as nib

    

def collect_genes(allexpressions_file, allcoord_file):
    """
    """

    human_home = "/usr/sci/projects/genes/human"
    # sub_dirs = ['donor9861', 'normalized_microarray_donor10021', 'normalized_microarray_donor14380', 'normalized_microarray_donor15697', 'normalized_microarray_donor12876', 'normalized_microarray_donor15496']
    sub_dirs = ['donor9861', 'normalized_microarray_donor10021']
    ontology_file_path = os.path.join(human_home, sub_dirs[0], 'Ontology.csv')
    probes_file_path = os.path.join(human_home, sub_dirs[0], 'Probes.csv')
    # mask_id = '4007'
    mask_id = '4005'

    # create a dictionary of ontology.
    gray_dic = {}
    with open(ontology_file_path, 'rb') as f:
        reader = csv.reader(f)
        reader.next() # skip header
        for row in reader:
            if mask_id in row[6]:
                gray_dic[row[0]] = True
            else:
                gray_dic[row[0]] = False
    
    expressions_all = list()
    coord_all = list()
    gene_filter = np.zeros(58692).astype(bool)
    for sub_dir in sub_dirs:
        print 'Working on {0}'.format(sub_dir)
        anno_file_path = os.path.join(human_home, sub_dir, 'SampleAnnot.csv')
        sample_filter = make_sample_filter(anno_file_path, gray_dic)
        print '    Sample filter size: {0} / {1}'.format(len([f for f in sample_filter if f == True]), len(sample_filter))
        coord_file = os.path.join(human_home, sub_dir, 'region_coordinates.txt')
        coord = np.loadtxt(coord_file)
        coord = coord[sample_filter, :]
        coord_all.append(coord)
        
        expression_file = os.path.join(human_home, sub_dir, 'MicroarrayExpression.csv')
        expressions = np.genfromtxt(expression_file, delimiter = ',')
        expressions = expressions[:,1:] # remove 1st column.
        expressions = expressions.T # now it's samples by genes.
        expressions = expressions[sample_filter,:]
        print '    expressions size after filtering samples: {0}'.format(expressions.shape)
        
        # compute gene filter.
        n_samples = expressions.shape[0]
        pacall_file = os.path.join(human_home, sub_dir, 'PACall.csv')
        gene_filter_sub = make_gene_filter(pacall_file, n_samples * 0.01, n_samples * 0.99)
        
        # keep a gene if it is good in any subject. We should be conservative.
        gene_filter = gene_filter | gene_filter_sub
        print '    gene_filter good genes: {0}'.format(len([i for i in gene_filter if i == True]))
        expressions_all.append(expressions)

    coord_all = np.vstack(coord_all)
    expressions_all = np.vstack(expressions_all)
    expressions_all = expressions_all[:, gene_filter] # filtering genes for all sub.
    print 'expressions size after filtering genes: {0}'.format(expressions_all.shape)
    np.save(allexpressions_file, expressions_all)
    np.savetxt(allcoord_file, coord_all, fmt = '%d')
    return coord_all, expressions_all

def make_sample_filter(sample_anno_file, dic):
    """
    make spatial sample filter.
    """

    filter = list()
    sample_coord = list()
    with open(sample_anno_file, 'rb') as f:
        reader = csv.reader(f)
        reader.next()
        for row in reader:
            filter.append(dic[row[0]])
            sample_coord.append([row[10], row[11], row[12]])

    return np.array(filter)

def make_gene_filter(pacall_file, th_low, th_high):
    """
    filter genes by expression.

    Some expressions are overly expressed, some others are not expressed at all. filter each row of the data matrix (genes x samples) so only moderately expressed remains.
    input: the PACall.csv, downloaed from the Allen institute website.
    output: a boolean np.ndarray that have True for included genes, and False for exclued.
    
    """
    
    pacall = np.genfromtxt(pacall_file, delimiter = ',', dtype = 'int')
    pacall = pacall[:,1:]
    rsum = pacall.sum(axis = 1)
    filter = (rsum < th_high) & (rsum > th_low)
    return np.array(filter)

def collect_fmri(coord_file, fmri_file):
    """
    collect fmri BOLD in all ROIs in one subject, and return a matrix of samples by time points.
    """

    coord = np.loadtxt(coord_file) # physical coordinates.
    fmri = nib.load(fmri_file)
    T = fmri.get_affine()
    S = np.linalg.inv(T[0:3, 0:3]) # rotation matrix.
    offset = T[0:3, 3] # shifting
    n_samples = coord.shape[0]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    D = fmri.get_data()
    n_timepoints = fmri.shape[3]
    fmri_mat = np.zeros((n_samples, n_timepoints))
    print 'fmri shape: {0}. Data type: {1}'.format(D.shape, fmri.get_data_dtype() )
    
    for wid in range(n_samples):
        v_idx = np.dot(S, coord[wid,:] - offset)
        v_idx = tuple(np.around(v_idx).astype(int))

        cand = [(i, j, k) for i in range(v_idx[0] - radius, v_idx[0] + radius + 1) for j in range(v_idx[1] - radius, v_idx[1] + radius + 1) for k in range(v_idx[2] - radius, v_idx[2] + radius + 1) if (i-v_idx[0])**2 + (j-v_idx[1])**2 + (k-v_idx[2])**2 <= radius**2]
        for v in cand:
            fmri_mat[wid,:] += D[v]

        fmri_mat[wid,:] = fmri_mat[wid,:] / len(cand)

    return fmri_mat

def collect_fmri_allsub(coord_file, session_dir, out_dir):
    """
    colletc fmri of all subjects, and save into a matrix for each sub.
    """
    sub_dirs = os.listdir(session_dir)
    sub_dirs.sort()
    for sub in sub_dirs:
        fmri_path = os.path.join(session_dir, sub, 'func/lfo_res2standard.nii.gz')
        fmri_mat = collect_fmri(coord_file, fmri_path)
        out_filename = os.path.join(out_dir, sub)
        np.save(out_filename, fmri_mat)
        print '{0} saved.'.format(out_filename)

def merge_subject_fmri(subject_dir, out_file):
    """
    merge npy file of all subjects. Run this routine after collect_fmri_allsub.
    """
    sub_files = os.listdir(subject_dir)
    sub_files.sort()
    bold_all = list()
    for sub in sub_files:
        bold_all.append(np.load(os.path.join(subject_dir, sub)))
        print 'sub {0} done.'.format(sub)

    bold_all = np.hstack(bold_all)
    np.save(out_file, bold_all)
        
    return bold_all

def test_pacall(pacall, mask_id, th_low, th_high):
    """
    test the truncating with different set of brain region.
    """
    human_home = "/usr/sci/projects/genes/human"
    ontology_file_path = os.path.join(human_home, 'donor9861', 'Ontology.csv')
    probes_file_path = os.path.join(human_home, 'donor9861', 'Probes.csv')

    # create a dictionary of ontology.
    gray_dic = {}
    with open(ontology_file_path, 'rb') as f:
        reader = csv.reader(f)
        reader.next() # skip header
        for row in reader:
            if mask_id in row[6]:
                gray_dic[row[0]] = True
            else:
                gray_dic[row[0]] = False

    gene_filter = np.zeros(58692).astype(bool)

    anno_file_path = os.path.join(human_home, 'donor9861', 'SampleAnnot.csv')
    sample_filter = make_sample_filter(anno_file_path, gray_dic)
    print '    Sample filter size: {0} / {1}'.format(len([f for f in sample_filter if f == True]), len(sample_filter))

    # compute gene filter.    
    n_samples = len([i for i in sample_filter if i == True])
    pacall = pacall[:, sample_filter]
    rmean = pacall.mean(axis = 1)
    plt.hist(rmean, 100)
    gene_filter = (rmean <= th_high) & (rmean >= th_low)
    
    print '    gene_filter good genes: {0}'.format(len([i for i in gene_filter if i == True]))

            
def collect_genes2(onto_id = '4005',
                    by_percent = False,
                    thr_exp_low = 0.01,
                    thr_exp_high = 0.99,
                    rm_noname_genes = True,
                    rm_XY = True,
                    log_xform = True,
                    ):
    """
    filtering genes and samples, after talk to Preeti.
    """
    human_home = "/usr/sci/projects/genes/human"
    sub_dirs = ['donor9861', 'normalized_microarray_donor10021', 'normalized_microarray_donor14380', 'normalized_microarray_donor15697', 'normalized_microarray_donor12876', 'normalized_microarray_donor15496']
    # sub_dirs = ['donor9861']
    ontology_file = "/usr/sci/projects/genes/human/donor9861/Ontology.csv"
    probes_file = os.path.join(human_home, sub_dirs[0], 'Probes.csv')
    n_genes = 58692

    # create a dictionary of ontology.
    gray_dic = {}
    with open(ontology_file, 'rb') as f:
        reader = csv.reader(f)
        reader.next() # skip header
        for row in reader:
            if onto_id in row[6]:
                gray_dic[row[0]] = True
            else:
                gray_dic[row[0]] = False

    expressions_all = list()
    coord_all = list()
    gene_filter = np.ones(n_genes).astype(bool)
    gene_filter_percent = np.zeros(n_genes).astype(bool)

    for sub_dir in sub_dirs:
        print 'Working on {0}'.format(sub_dir)
        anno_file_path = os.path.join(human_home, sub_dir, 'SampleAnnot.csv')
        sample_filter = make_sample_filter(anno_file_path, gray_dic)
        print '    Sample filter size: {0} / {1}'.format(len([f for f in sample_filter if f == True]), len(sample_filter))
        coord_file = os.path.join(human_home, sub_dir, 'region_coordinates.txt')
        coord = np.loadtxt(coord_file)
        coord = coord[sample_filter, :]
        coord_all.append(coord)

        # first check if there is .npy file.
        expression_npy = os.path.join(human_home, sub_dir, 'preprocessing', 'MicroarrayExpression.npy')
        if os.path.exists(expression_npy):
            expressions = np.load(expression_npy)
            print '    expressions file loaded from {0}.'.format(expression_npy)
        else:
            expression_file = os.path.join(human_home, sub_dir, 'MicroarrayExpression.csv')
            expressions = np.genfromtxt(expression_file, delimiter = ',')            

        expressions = expressions[:,1:] # remove 1st column.
        expressions = expressions.T # now it's samples by genes.
        expressions = expressions[sample_filter,:]
        print '    expressions size after filtering samples: {0}'.format(expressions.shape)
        
        if by_percent:
            n_samples = expressions.shape[0]
            pacall_npy = os.path.join(human_home, sub_dir, 'preprocessing', 'PACall.npy')
            if os.path.exists(pacall_npy):
                pacall = np.load(pacall_npy)

            else:
                pacall_file = os.path.join(human_home, sub_dir, 'PACall.csv')
                pacall = np.genfromtxt(pacall_file, delimiter = ',', dtype = 'int')
                
            pacall = pacall[:,1:]                
            rsum = pacall.sum(axis = 1)
            gene_filter_percent_sub = (rsum < n_samples * thr_exp_high) & (rsum > n_samples * thr_exp_low)
            gene_filter_percent = gene_filter_percent | gene_filter_percent_sub
            print '    gene_filter_percent good genes: {0}'.format(len([i for i in gene_filter_percent if i == True]))

        expressions_all.append(expressions)
        
    coord_all = np.vstack(coord_all)
    expressions_all = np.vstack(expressions_all)

    # remove over/under expressed ?
    if by_percent:
        gene_filter = gene_filter & gene_filter_percent

    # remove useless (noname) probes?
    if rm_noname_genes:
        gene_filter_useful = list()
        with open(probes_file, 'rb') as f:
            reader = csv.reader(f)
            reader.next() # skip header
            for row in reader:
                if row[5]:
                    # entrez_id not empty
                    gene_filter_useful.append(True)
                else:
                    gene_filter_useful.append(False)
                    
        gene_filter_useful = np.asarray(gene_filter_useful)
        gene_filter = gene_filter & gene_filter_useful
        print 'After removing useless,  good genes: {0}'.format(len([i for i in gene_filter if i == True]))

    # remove probles of X or Y genes?
    if rm_XY:
        gene_filter_XY = list()
        with open(probes_file, 'rb') as f:
            reader = csv.reader(f)
            reader.next()
            for row in reader:
                if (not "X" in row[6]) & (not "Y" in row[6]) :
                    gene_filter_XY.append(True)
                else:
                    gene_filter_XY.append(False)

        gene_filter_XY = np.asarray(gene_filter_XY)
        gene_filter = gene_filter & gene_filter_XY
        print 'After removing XY chromosome,  good genes: {0}'.format(len([i for i in gene_filter if i == True]))

    # finally apply the gene filter.
    expressions_all = expressions_all[:, gene_filter] # filtering genes for all sub.

    # do log transformation
    if log_xform:
        expressions_all = np.log(expressions_all)
        
    return coord_all, expressions_all, gene_filter
            
