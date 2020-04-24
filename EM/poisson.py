#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 01:57:11 2020

@author: jaeykim
EM for mixed model of (pois(lambda) with drop off chance of alpha)
"""

from scipy.stats import poisson
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import fsolve
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from scipy.special import logsumexp
import numpy as np
import math
import os
import matplotlib.pyplot as plt


"""
# solves for x in e^x(u - ax) = u

example
sol = print(numericSolveLinearExp(10, 5, tol=1e-5))
Note that u >= a
    
"""

def print_classification(true_classification, pred, num_K):
    min_perf = 1
    total_perf = 0
    for i in range(num_K):
        counts = np.bincount(true_classification[pred == i])
        perf = np.max(counts) / sum(pred == i)
        #print("cluster", i,"correct:", perf)
        if perf < min_perf:
            min_perf = perf
        total_perf += perf
    print("ARI:", adjusted_rand_score(true_classification, pred),
          "min perf:", min_perf, "avg perf:", total_perf / num_K)
    
def numericSolveLinearExp(u, a, tol = 1e-4):
    high = u / a
    low = 0
    #assert(high + np.log(u - a * high) < np.log(u))
    #assert(low + np.log(u - a * low) >= np.log(u))
    
    guess = (low + high) / 2
    while True: 
        if guess + np.log(u - a * guess) > np.log(u): 
            low = guess 
        else: 
            high = guess 
        guess = ( low + high ) / 2
        if high - low < tol:
            break
    
    return guess

# return dropoff, mu (mean)
def MLE(n, k, u):
    if k == n:
        assert(u == 0)
        return 1, 0
    
    mu = u / (n - k)
    prev_mu = 0
    prev_dropoff = k / n
    dropoff = k / n
    while abs(prev_mu - mu) > 1e-6 or abs(dropoff - prev_dropoff) > 1e-6:
        prev_mu = mu
        prev_dropoff = dropoff
        mu = u / (n - n * dropoff)
        dropoff = max(k - n * np.exp(-mu), 0) / n # first guess
        
    assert(dropoff != 1)
        
    return dropoff, mu

def init_cluster(X, num_K, num_neighbors):
    print("clustering")
    [num_rows, num_cols] = np.shape(X)
    X_logged = np.log(X + 1)
    X_normalized = X_logged / np.std(X_logged, axis = 1).reshape(num_rows, 1)
    clustering = SpectralClustering(n_clusters=num_K, n_neighbors = num_neighbors, affinity='nearest_neighbors')
    label_spectral = clustering.fit_predict(X_normalized)
    print("finished clustering")
    return label_spectral

def normalize_data(A):
    #  Simple convenience function to library and log normalize a matrix
    # column sum
    totalUMIPerCell = A.sum(axis=0)

    A_norm = A / totalUMIPerCell
    A_norm = A_norm * 1e4
    A_norm = np.log(A_norm +1)
    return A_norm
    
def kmeans_cluster(X, num_K):
    X_norm = normalize_data(X.T).T
    X_pca = PCA(n_components = 100).fit_transform(X_norm)
    labels = KMeans(n_clusters = num_K).fit_predict(X_pca)
    return labels


def MixedPoisson(X, num_K, num_iter = 10, true_labels = []):
    """
        input trainX is a N by D matrix containing N datapoints, num_K is the 
        number of clusters or mixture components desired.
        num_iter is the maximum number of EM iterations run over the dataset
        Description of other variables:
            - mu which is K*D, the coordinates of the means
            - pk, which is K*1 and represents the cluster proportions
            - zk, which is N*K, has at each z(n,k) the probability that the 
            nth data point belongs to cluster k, specifying the cluster 
            associated with each data point
            - si2 is the estimated (shared) variance of the data
            - BIC is the Bayesian Information Criterion 
            (smaller BIC is better)
    """
    N = X.shape[0]
    D = X.shape[1]

    try:
        if num_K >= N:
            raise AssertionError
    except AssertionError:
        print("You are trying too many clusters")
        raise


    
    pk = np.ones((num_K,1))/num_K # Uniformly initialize cluster proportions
    dropoffs = np.zeros([D,num_K])
    mu = np.ones([D, num_K]) # Random initialization of clusters    
    
    # use hard assignment to reduce numeric instabilities
    #labels = true_labels - 1
    #num_true = 500
    #labels = np.random.randint(0, num_K, N)
    #labels[:num_true] = true_labels[:num_true] - 1
    #labels = init_cluster(X, num_K, 5) 
    labels = kmeans_cluster(X, num_K)
            
    if len(true_labels) != 0:
        print("Initial clustering")
        print_classification(true_labels, labels, num_K)
    
    
    def E_step():
        loglike = 0
        oldlike = 0
        num_labels_changed = 0
        for cell in range(N):
            zero_idx = X[cell,:] == 0
            non_zero = X[cell,:] > 0
            logprob1 = [np.sum(np.log(dropoffs[zero_idx, k] + (1 - dropoffs[zero_idx, k]) * np.exp(-mu[zero_idx, k]))) for k in range(num_K)]
            logprob2 = [np.sum(np.log(1 - dropoffs[non_zero, k])) for k in range(num_K)]
            logprob3 = [np.sum(poisson(mu[non_zero,k]).logpmf(X[cell,non_zero])) for k in range(num_K)]                   
            
            logprob = np.array(logprob1) + np.array(logprob2) + np.array(logprob3)
            assert(not (logprob == float("-inf")).all())
            new_label = np.argmax(logprob)
            loglike += logprob[new_label]
            oldlike += logprob[labels[cell]]
            if new_label != labels[cell]:
                num_labels_changed = num_labels_changed + 1
                labels[cell] = new_label
        print("Changed", num_labels_changed, "labels")
        print("new log like", loglike, "old", oldlike)
        
        return num_labels_changed
            
    
    def M_step():        
        dropoffs_tol = 0
        mu_lower = 0
        # optimal value of lambda in cluster e^lambda(U - (n - k) * lambda) = u
        # U sum of members, n num members, k num zeros
        for k in range(num_K):
            num_members = sum(labels == k)
            pk[k] = num_members / N
            for gene in range(D):
                num_zeros = sum(X[labels == k, gene] == 0)
                sum_members = sum(X[labels == k, gene])
                dropoffs[gene, k], mu[gene, k] = MLE(num_members, num_zeros, sum_members)

        dropoffs[dropoffs < dropoffs_tol] = dropoffs_tol
        dropoffs[dropoffs > 1 - dropoffs_tol] = 1 - dropoffs_tol
        mu[mu < mu_lower] = mu_lower
        """
        for cell in range(N):
            assert(not(dropoffs[X[cell, :] > 0, labels[cell]] == 1).any())
        """
        
    for i in range(num_iter):
        print("starting iteration", i)
        M_step()
        num_changed = E_step()
        print("pk", pk)
        #print("mu", mu)
        print("dropoffs 0 or 1", np.sum(dropoffs == 0) + np.sum(dropoffs == 1))
        if len(true_labels) != 0:
            print_classification(true_labels, labels, num_K)
        if num_changed == 0:
            break
        
        
    print("mu", mu)
    print("dropoffs", dropoffs)
    pred = X.copy()
    pred = pred.astype(float)
    for cell in range(N):
        pred[cell, X[cell,:] == 0] = mu[X[cell,:] == 0, labels[cell]]

    return labels, pred
        

def load_helper(subfolder, fname):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    try:
        file_to_load = data_folder + fname
        with open(file_to_load) as data_file:
#            ncols = len(data_file.readline().split(','))
#            X = np.genfromtxt(data_file,dtype = 'int_',delimiter=',',usecols=range(1,ncols))
            X = np.genfromtxt(data_file,dtype = 'int_',delimiter=' ') #TODO change back to int
        Y = np.matrix(X)
    except FileNotFoundError:
        file_to_load = data_folder + fname
        with open(file_to_load) as data_file:
#            ncols = len(data_file.readline().split(','))
#            X = np.genfromtxt(data_file,delimiter=',',usecols=range(1,ncols))
            X = np.genfromtxt(data_file,delimiter=' ')
        Y = np.matrix(X)
    return Y

def load_counts(subfolder):
    return load_helper(subfolder, 'counts.csv')

def load_true_counts(subfolder):
    return load_helper(subfolder, 'true_counts.csv')

def load_classification(subfolder):
    return load_helper(subfolder, 'group_data.csv')

# [input]: subfolder of SplatGenData where csv lives, name of csv. Also can load .csv.gz (g-zipped csv)
# behaviour: skips the first row and first column
def load_csv_to_matrix(subfolder,name):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    file_name = data_folder + name
    with open(file_name) as data_file:
#        ncols = len(data_file.readline().split(','))
#        X = np.genfromtxt(data_file,dtype = 'int_',delimiter=',',usecols=range(1,ncols))
        X = np.genfromtxt(data_file,dtype = 'int_',delimiter=' ')
    Y = np.matrix(X)
    return Y
    
def MSE(X_true, X_pred, X_observed):
    return np.mean((np.array(X_true - X_pred) ** 2))

def MSE_restricted(X_true, X_pred, X_observed):
    return np.mean((np.array(X_true[X_observed == 0] - X_pred[X_observed == 0]) ** 2))


    
"""
Here row is cells column are genes
"""
def main():
    """
    n = 1000
    dropoff = 0.8
    mu = 3
    k = (dropoff + np.exp(-mu)) * n
    u = n * (1 - dropoff) * mu
    MLE(n,k,u)
    return 0
    """
    
    
    num_K = 10
    print("loading data")
    #subfolder = 'new_small_examples/'
    subfolder = '10_groups_1000_cells_5000_genes/';
    X = np.array(load_counts(subfolder).T)
    X_true = np.array(load_true_counts(subfolder).T)
    [num_cells, num_genes] = np.shape(X)
    print("finished loading")
    
    true_labels = np.array(load_classification(subfolder)).reshape(num_cells)

    labels, pred = MixedPoisson(X, num_K, 1, true_labels = true_labels)
    print("MSE_restricted", MSE_restricted(X_true, pred, X))
    print(sum(pred > 100))
    #print(((X_true - pred)[X == 0])[0:500])
    #print(((X_true - X)[X == 0])[0:500])


if __name__ == '__main__':
    main()

    
            
            
        
                
