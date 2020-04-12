#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:44:17 2020

@author: jaeykim
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import os

class VecCluster:
    # n number of rows
    def __init__(self, n):
        self.n = n # number of rows
        # number of vectors that is 1 for this entry
        # sv = sum of vectors added
        self.pop = np.zeros(n)
        # number of vectors thats is assigned and covered by this vec
        self.num_vec = 0
        
    def SV(self):
        return 1 * (self.pop > 0)
        
    def score(self):
        return self.num_vec * sum(self.pop > 0) - sum(self.pop)
        
    # Expects a 0-1 vector in  {0,1}^n
    # Returns the increase in score
    def addVec(self, V):
        old_score = self.score()
        self.pop += V
        self.num_vec += 1
        return self.score() - old_score
        
    # One can only remove vectors that was added before!
    def removeVec(self, V):
        self.pop -= V
        self.num_vec -= 1
        self.checkInvariant()
        
    # Returns the increase in score
    def addVecCluster(self, VC):
        old_score = self.score()
        self.pop += VC.pop
        self.num_vec += VC.num_vec
        return self.score() - old_score
        
    def removeVecCluster(self, VC):
        self.pop -= VC.pop
        self.num_vec -= VC.num_vec
        self.checkInvariant()
        
    def checkInvariant(self):
        assert(self.num_vec >= 0)
        assert((self.pop >= 0).all())

def improve_clustering(X, rank, init_label = [], 
                               num_steps = float("inf")):
    [r,c] = np.shape(X)
    SupportVectors = [VecCluster(r) for i in range(rank)]
    def score():
        return sum([sv.score() for sv in SupportVectors])
    
    labels = []
    if len(init_label) == 0:
        labels = [-1 for i in range(c)]
    elif len(init_label) == c:
        labels = init_label
        for i in range(c):
            SupportVectors[labels[i]].addVec(X[:,i])
        print("initial score", score())
    else:
        print("ERROR! expected init_label size of",c, 
              "but got", len(init_label))
        
    
    #shift vectors
    #returns number of vectors moved
    def shift_vec():
        num_changed = 0
        for i in range(c):
            if labels[i] != -1:
                SupportVectors[labels[i]].removeVec(X[:,i])
                
            min_increase_in_score = float("inf")
            min_idx = labels[i] # this makes it so that we don't change labels often
            for sv in range(rank):
                increase_in_score = SupportVectors[sv].addVec(X[:,i])
                SupportVectors[sv].removeVec(X[:,i])
                if increase_in_score < min_increase_in_score:
                    min_increase_in_score = increase_in_score
                    min_idx = sv
            SupportVectors[min_idx].addVec(X[:,i])
            if min_idx != labels[i]:
                labels[i] = min_idx
                num_changed += 1
                
        return num_changed
            
    #shift sub clusters (all vec in a cluster that has a 1 on a certain row)
    #returns number of sub clusters moved
    def shift_sub_cluster():
        num_changed = 0
        for i in range(rank):
            for row in range(r):
                sub_vec_cluster = VecCluster(r)
                sub_vec_cluster_idx = []
                for j in range(c):
                    if X[row,j] == 1 and labels[j] == i:
                        sub_vec_cluster.addVec(X[:,j])
                        sub_vec_cluster_idx.append(j)
                if len(sub_vec_cluster_idx) == 0:
                    continue
                
                SupportVectors[i].removeVecCluster(sub_vec_cluster)
                        
                min_increase_in_score = float("inf")
                min_cluster = i
                for sv in range(rank):
                    increase_in_score = SupportVectors[sv].addVecCluster(sub_vec_cluster)
                    SupportVectors[sv].removeVecCluster(sub_vec_cluster)
                    if increase_in_score < min_increase_in_score:
                        min_increase_in_score = increase_in_score
                        min_cluster = sv
                SupportVectors[min_cluster].addVecCluster(sub_vec_cluster)
                if min_cluster != i:
                    # We need to move!
                    for idx in sub_vec_cluster_idx:
                        labels[idx] = min_cluster
                    num_changed += 1
        return num_changed
            
    vec_stopped = False
    cluster_stopped = False
    num_iter_finished = 0
    while True:
        vec_changed = shift_vec()
        print("shifted", vec_changed,"vectors score:",score())
        vec_stopped = vec_changed == 0
        if vec_stopped and cluster_stopped:
            break
        cluster_changed = shift_sub_cluster()
        print("shifted", cluster_changed,"sub clusters score:", score())
        cluster_stopped = cluster_changed == 0
        if vec_stopped and cluster_stopped:
            break
        num_iter_finished += 1
        if num_iter_finished == num_steps:
            print("Terminating early due to reaching maximum number of iterations")
        
    return SupportVectors, labels

def test_cluster_span_deterministic():
    X = np.zeros([4,6])
    X[:,0] = 1
    X[:, 1] = [1,1,0,0]
    X[:, 2] = [0,0,1,1]
    X[:,3] = [0,1,1,0]
    X[:,4] = X[:,3]
    X[:,5] = X[:,3]
    improve_clustering(X, 2)
    

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

''' 
WARNING: As of now, your working directory should be the TechnicalZeroDetector/ 
folder in order for the load function to work correctly.
'''

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

def load_counts(subfolder):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    try:
        file_to_load = data_folder + 'counts.csv'
        with open(file_to_load) as data_file:
            ncols = len(data_file.readline().split(','))
            X = np.genfromtxt(data_file,dtype = 'int_',delimiter=',',usecols=range(1,ncols))
        Y = np.matrix(X)
    except FileNotFoundError:
        file_to_load = data_folder + ' counts.csv'
        with open(file_to_load) as data_file:
            ncols = len(data_file.readline().split(','))
            X = np.genfromtxt(data_file,delimiter=',',usecols=range(1,ncols))
        Y = np.matrix(X)
    return Y

# [input]: subfolder of SplatGenData where csv lives, name of csv. Also can load .csv.gz (g-zipped csv)
# behaviour: skips the first row and first column
def load_csv_to_matrix(subfolder,name):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    file_name = data_folder + name
    with open(file_name) as data_file:
        ncols = len(data_file.readline().split(','))
        X = np.genfromtxt(data_file,dtype = 'int_',delimiter=',',usecols=range(1,ncols))
    Y = np.matrix(X)
    return Y

# [input]: classification vector and subfolder of SplatGenData to which it belongs
# [output]: 1 if successful?
#writes classification as classification.csv in the desired subfolder
def write_classification(classification,subfolder):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    file_name = data_folder + 'classification.csv'

# [input]: mask, subfolder of SplatGenData to which it belongs, name to save (usually "tight" or "loose") including extension. Allows gz zipping (e.g. give name "tightmask.csv.gz").
# [output]: 1 if successful?
def write_matrix_to_data(matrix,subfolder,name):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    file_name = data_folder + name
    np.savetxt(file_name,matrix,delimiter=',')


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
    TightMask, LooseMask, classification = cluster_span(X_corrupted.copy(), rank, 0.07)

    SV, SVclassification = improve_clustering(X_corrupted.copy(), 
                                                      rank, classification.astype(int))
    
    print("finished deter")
    X_SV = np.zeros([row, col])
    for i in range(col):
        X_SV[:,i] = SV[SVclassification[i]].SV()
        
    
    
    fig, axs = plt.subplots(5)
    fig.suptitle('0-1 reconstruction')
    axs[0].imshow(X)
    axs[1].imshow(X_corrupted)
    axs[2].imshow(LooseMask)
    axs[3].imshow(TightMask)
    axs[4].imshow(X_SV)
    plt.savefig("0-1_reconstruction", dpi=600)


def main():
    subfolder = '5_groups_1000_cells_5000_genes/';
    X = load_counts(subfolder)
    [row,col] = np.shape(X)
    
    # Binarize counts matrix
    X_binary = np.where(X>0, 1,0)

    threshold_ratio = 0.2  # TODO take as argument see description in documentation
    rank = 5 # TODO take as an argument (should be an upper bound on clustering)

    # Run the technical zero finding algorithm
    TightMask, LooseMask, classification = cluster_span(X_binary, rank, threshold_ratio)

    # Save TightMask, LooseMask, Classification
    write_matrix_to_data(TightMask,subfolder,'tight_mask.csv.gz')
    write_matrix_to_data(LooseMask,subfolder,'loose_mask.csv.gz')
    write_matrix_to_data(classification,subfolder,'classification.csv')

    # TODO print some performance stuff? Number of technical zeros detected (hopefully)

    # \Jeremy: I'm not going to write the matrix with this reordering because I want the
    # saved matrix to have the same cell-ordering as the input. But we can make figures now
    # and they'll look reordered by classification
    sorted_indices = np.argsort(classification)
    idx = np.empty_like(sorted_indices)
    idx[sorted_indices] = np.arange(len(sorted_indices))
    TightMask[:] = TightMask[:, idx]
    LooseMask[:] =LooseMask[:, idx]
    classification[:] = classification[sorted_indices]
    
    # TODO: Should try to compute differentially expressed genes (genes that differ btwn cell types)
    # and only plot those .
    # TODO: load ground truth (if available) and also display that.
    fig, axs = plt.subplots(4)
    fig.suptitle('0-1 reconstruction')
    axs[0].imshow(classification.reshape(col,1).T)
    axs[1].imshow(X_binary)
    axs[2].imshow(LooseMask)
    axs[3].imshow(TightMask)
    plt.savefig("0-1_reconstruction_jp", dpi=600)


if __name__ == '__main__':
    main()
