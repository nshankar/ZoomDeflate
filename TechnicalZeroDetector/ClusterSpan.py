#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:44:17 2020

@author: jaeykim
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

"""
[Input]
X: 0-1 matrix. Column for each cell, row for each gene
    X[g,c] = 1 indicates that gene g was observed in cell c
rank: Upper bound on number of cell types/clusters expected
threshold_ratio: this is a hyper-parameter. Ensure that the threshold_ratio
    is less than the 1 - dropoff rate. 
    E.g. for sparcity 50% following seemed to work well
    you need to play around with it
    dropoff 50% threshold_ratio = 0.4
    dropoff 70% threshold_ratio = 0.2
    
[Output]
[TightMask, LooseMask, classification]
Masks - 0 indicates that we predict that it is a true zero
TightMask is more likely to output 0 
    (less chance of falsely classifying techinical zero as a biological zero)
LooseMask is more likely to output 1 
    (more likely to detect technical zeros)

cluster_span finds rank many vectors V_1, V_2,..., V_rank 
that spans the image of the matrix X
It ensures that for each cell, there is a V_i such that genes observed
in the cell is a subset of genes corresponding to V_i.
This function attempts to minimizes the number of zeros of each spanning vectors
such that each cell fits in very "tight" to some spanning vector V_i


"""
def cluster_span(X, rank, threshold_ratio):
    [r,c] = np.shape(X)
    # Run kmeans to get a general idea of clustering
    kmeans = KMeans(n_clusters=rank).fit(X.T)
    
    # Get some genes we are confident that belongs to the clusters
    # this partially initializes spanning vectors
    clusters = [X[:,kmeans.labels_ == i] for i in range(rank)]
    span_vec = [np.zeros(r) for i in range(rank)]
    non_one_val = -1
    for i in range(rank):
        [r_c, c_c] = np.shape(clusters[i])
        threshold = threshold_ratio * c_c
        for j in range(r):
            span_vec[i][j] = 1 if sum(clusters[i][j,:]) >= threshold else non_one_val
                
    # Now complete the spanning vectors
    X[X == 0] = non_one_val
    for i in range(c):
        best_vec = -1
        best_score = float("-inf")
        # Try to fit it in to a span_vec first
        for j in range(rank):
            if not (X[:,i] <= span_vec[j]).all():
                continue
            
            score = np.dot(X[:,i], span_vec[j])
            if score > best_score:
                best_score = score
                best_vec = j
            
        if best_vec == -1:
            for j in range(rank):
                score = np.dot(X[:,i], span_vec[j])
                if score > best_score:
                    best_score = score
                    best_vec = j
            
        # Now "merge" cell i into to best span_vec
        for k in range(r):
            if X[k,i] == 1 and span_vec[best_vec][k] == non_one_val:
                span_vec[best_vec][k] = 1
                
    # Now classify points using spanning vectors and complete X
    LooseMask = np.zeros([r,c])
    TightMask = np.ones([r,c])
    
    classification = np.zeros(c)
    for i in range(c):
        #TODO below for loop is duplicated. Can I avoid passing by value?
        best_vec = -1
        best_score = float("-inf")
        for j in range(rank):
            if not (X[:,i] <= span_vec[j]).all():
                continue
            
            TightMask[span_vec[j] == non_one_val,i] = 0
            
            score = np.dot(X[:,i], span_vec[j])
            if score > best_score:
                best_score = score
                best_vec = j
        if best_vec == -1:
            print("ERROR! Spanning set does not actually span")
        classification[i] = best_vec
        LooseMask[:,i] = span_vec[best_vec]
        
    return TightMask, LooseMask, classification
        

#random 0-1 mat of rank at most rank
#Warning! not very random, just for testing
def random_mat(r,c,rank,sparcity):
    n = int(np.floor(c/rank))
    X = np.zeros([r,c])
    for i in range(rank):
        rand_col = 1 * (np.random.rand(r,1) > sparcity)
        for j in range(n):
            X[:,i * n + j] = rand_col.T
            
    return X

def test():
    row = 50
    col = 450
    rank = 6
    dropoff = 0.7
    sparcity = 0.5
    X = random_mat(row, col, rank, sparcity)
    cor_mat = np.random.rand(row, col)
    X_corrupted = np.zeros([row, col])
    X_corrupted[cor_mat > dropoff] = X[cor_mat > dropoff]
    TightMask, LooseMask, classification = cluster_span(X_corrupted, rank, 0.2)
    
    kmeans = KMeans(n_clusters=rank).fit(X_corrupted.T)
    X_weighted = X_corrupted.copy()
    X_weighted[X_weighted == 0] = -0.5
    
    kmeans_weighted = KMeans(n_clusters=rank).fit(X_weighted.T)
    
    
    print(np.shape(TightMask))
    
    fig, axs = plt.subplots(6)
    fig.suptitle('0-1 reconstruction')
    axs[0].imshow(kmeans.labels_.reshape(col,1).T)
    axs[1].imshow(classification.reshape(col,1).T)
    axs[2].imshow(X)
    axs[3].imshow(X_corrupted)
    axs[4].imshow(LooseMask)
    axs[5].imshow(TightMask)
    plt.savefig("0-1_reconstruction", dpi=600)

def main():
    X = # TODO import 0-1 matrix
    threshold_ratio = # TODO take as argument see description in documentation
    rank = # TODO take as an argument (should be an upper bound on clustering)
    
    TightMask, LooseMask, classification = cluster_span(X, rank, threshold_ratio)
    
    # TODO save TightMask LooseMask
    # TODO print some performance stuff? Number of technical zeros detected (hopefully)
    # Also, it might be nice to see what happens when we reorder X and masks
    # based on classification of cells to see how things line up
    
    

if __name__ == '__main__':
    main()