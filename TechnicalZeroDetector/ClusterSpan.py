#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:44:17 2020

@author: jaeykim
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def reconstruct(X, row_k, col_k,threshold = 0.5):
    col_kmeans = KMeans(n_clusters=row_k).fit(X.T)
    row_kmeans = KMeans(n_clusters=col_k).fit(X)
    for row_label in range(row_k):
        for col_label in range(col_k):
            in_col = col_kmeans.labels_ == col_label
            in_row = row_kmeans.labels_ == row_label
            region = np.matmul(in_row.reshape(20,1), in_col.reshape(1,25))
            num_in_region = sum(sum(region))
            num_one = sum(X[region] == 1)
            if num_one >= threshold * num_in_region:
                X[region] = 1
                
    print("cols", col_kmeans.labels_)
    print("rows", row_kmeans.labels_)
    return X

# Gene ranking by average value of inner product of columns given
# it has this gene divided by avg inner prod when its off
def conditional_info_score(X):
    # TODO normalize!
    InnerProd = np.matmul(X.T, X)
    [r, c] = np.shape(X)
    score = np.zeros(r)
    for gene in range(r):
        has_gene = X[gene,:] != 0
        if sum(has_gene) <= 1:
            continue #set score to zero
        sim_gene_on = 0
        sim_gene_on_pair_count = 0
        sim_gene_off = 0
        sim_gene_off_pair_count = 0
        for i in range(c):
            for j in range(i+1,c):
                if not (has_gene[i] or has_gene[j]):
                    sim_gene_off += InnerProd[i,j]
                    sim_gene_off_pair_count += 1
                    #print("off",gene,i,j)
                elif has_gene[i] and has_gene[j]:
                    sim_gene_on += InnerProd[i,j]
                    sim_gene_on_pair_count += 1
                    #print("on",gene,i,j)
        score[gene] = (sim_gene_on /sim_gene_on_pair_count)\
            / (sim_gene_off / sim_gene_off_pair_count)
        
    return score
    
    
# k number of weights
def AdaBoost(X, col_label, k):
    
    [r, c] = np.shape(X)
    weights = np.ones(c) / c
    row_score = np.zeros(r)
    for i in range(r):
        row = X[i,:]
        adjusted_label = col_label
        adjusted_label[row == 0] *= -1
        one_weights = [sum(weights.T[adjusted_label == l]) 
                       for l in range(k)]
        # one_weights[l] is the sum of weights that are 1 on the row
        # and labeled l
        max_one_label = np.argmax(one_weights)
        zero_weights = [sum(weights[adjusted_label == -l]) 
                       for l in range(k)]
        max_zero_label = np.argmax(zero_weights)
        
        ## Compute error
        # Weights add up to 1 so just subtract weights we get correct
        err = 1 - one_weights[max_one_label] - zero_weights[max_zero_label]
        row_score[i] = np.log((1 - err) / err) + np.log(k - 1)
        
        ## Update Weights
        for j in range(c):
            # If we get the correct label, we don't need to update
            if row[j] == 0:
                if col_label[j] == max_zero_label:
                    continue
            else:
                if col_label[j] == max_one_label:
                    continue
                
            weights[j] = weights[j] * np.exp(row_score[i])
            
        ## Normalize weights
        weights = weights / sum(weights)
        
    print(row_score)
    return row_score
        
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

def min_span(X, rank):
    threshold_ratio = 0.4 # TODO fix this
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
    X_completed = np.zeros([r,c])
    classification = np.zeros(c)
    for i in range(c):
        #TODO below for loop is duplicated. Refactor!
        best_vec = -1
        best_score = float("-inf")
        #num_span_fit = 0
        for j in range(rank):
            if not (X[:,i] <= span_vec[j]).all():
                continue
            #num_span_fit += 1
            
            score = np.dot(X[:,i], span_vec[j])
            if score > best_score:
                best_score = score
                best_vec = j
        #print(num_span_fit)
        if best_vec == -1:
            print("error!")
        classification[i] = best_vec
        X_completed[:,i] = span_vec[best_vec]
        
    return classification, X_completed
        

def main():
    row = 50
    col = 450
    rank = 6
    dropoff = 0.5
    sparcity = 0.5
    X = random_mat(row, col, rank, sparcity)
    cor_mat = np.random.rand(row, col)
    X_corrupted = np.zeros([row, col])
    X_corrupted[cor_mat > dropoff] = X[cor_mat > dropoff]
    classification, X_completed = min_span(X_corrupted, rank)
    
    kmeans = KMeans(n_clusters=rank).fit(X_corrupted.T)
    X_weighted = X_corrupted.copy()
    X_weighted[X_weighted == 0] = -0.5
    
    kmeans_weighted = KMeans(n_clusters=rank).fit(X_weighted.T)
    
    
    
    
    fig, axs = plt.subplots(5)
    fig.suptitle('0-1 reconstruction')
    axs[0].imshow(kmeans.labels_.reshape(col,1).T)
    axs[1].imshow(classification.reshape(col,1).T)
    axs[2].imshow(X)
    axs[3].imshow(X_corrupted)
    axs[4].imshow(X_completed)
    plt.savefig("0-1_reconstruction", dpi=600)
    

if __name__ == '__main__':
    main()