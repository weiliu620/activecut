import sys
import os
import numpy as np
import pylab as pl
from time import time
import mdp

from sklearn import decomposition
from sklearn import  metrics
from sklearn.cluster import KMeans
from sklearn import preprocessing

def pca(X, n_comp):
    """
    regular pca.
    """
    estimator = decomposition.PCA(n_components = n_comp)
    estimator.fit(X)
    X_proj = estimator.transform(X)
    return estimator.components_, X_proj, estimator.explained_variance_ratio_

def skpca(X, n_comp):
    """
    PCA by scikit-learn.
    """
    pca = decomposition.PCA()
    pca.fit(X)

def kernelpca(X, n_comp):
    """
    kernel pca. But way too slow.
    """
    estimator = decomposition.KernelPCA(n_components = n_comp, kernel = 'rbf')
    estimator.fit(X)
    X_proj = estimator.transform(X)
    return estimator.components_, X_proj, 
    

def sparsepca(X, n_comp):
    """
    run sprase pca on data.
    """
    n_samples, n_features = X.shape
    # center the data. Note we only do the global centering. i.e. earch column
    # is centered to zero mean. Though the data should already 'locally'
    # centered since each row is normalized to z score. 
    X = X - X.mean(axis = 0)
    estimator = decomposition.SparsePCA(n_components=n_comp, alpha=0.8, max_iter = 100, n_jobs = 20, verbose = 1, tol = 1e-2)
    t0 = time()
    estimator.fit(X)
    train_time = (time() - t0)
    print "done in %0.3fs" % train_time
    components_ = estimator.components_
    X_projected = estormator.transform(X)
    
    return components_, X_projected

def sparsepca_minibatch(X, n_comp):
    """
    run sprase pca on data.
    """
    n_samples, n_features = X.shape
    estimator = decomposition.MiniBatchSparsePCA(n_components=n_comp, alpha=0.8, n_iter = 100, n_jobs = 2, verbose = True, batch_size = 3)
    t0 = time()
    estimator.fit(X)
    train_time = (time() - t0)
    print "done in %0.3fs" % train_time
    components_ = estimator.components_
    X_projected = estimator.transform(X)
    
    return components_, X_projected

def sparsepca_minibatch2(X, n_comp):
    """
    run sparse pca mini batch with MDP.
    """
        
    node = mdp.nodes.MiniBatchSparsePCAScikitsLearnNode(n_components = n_comp, alpha = 0.8,n_iter = 1000, batch_size = 3, verbose = 2)
    node.train(X)
    node.stop_training()
    print 'output dimension: {0}'.format(node.output_dim)
    X_projected = node.execute(X)
    return node.scikits_alg.components_, X_projected
                                         
def kmeans(X, n_clust):
    """
    kmeans clustering on X (samples x features)

    Note this is differnt from pca. The spatial locations are samples, and genes are features. Use expressions.T() as X.
    """

    X = scale(X)
    estimator = KMeans(init = 'k-means++', n_clusters = n_clust, n_init = 10, verbose = 2)
    
    estimator.fit(X)
    labels = estimator.predict(X)
    return labels
        
        
