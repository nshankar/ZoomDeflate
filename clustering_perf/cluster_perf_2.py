#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 19:15:10 2020

@author: jaeykim
"""

import numpy as np
from sklearn.cluster import KMeans
import os
from sklearn.decomposition import PCA
from numpy import percentile
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import PCA

def load_helper(subfolder, fname,header_lines):
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    data_folder = data_path + subfolder
    try:
        file_to_load = data_folder + fname
        with open(file_to_load) as data_file:
#            ncols = len(data_file.readline().split(','))
#            X = np.genfromtxt(data_file,dtype = 'int_',delimiter=',',usecols=range(1,ncols))
            X = np.genfromtxt(data_file,dtype = 'int_',delimiter=' ',skip_header=header_lines) #TODO change back to int
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
    return load_helper(subfolder, 'counts.csv',0)

def load_true_counts(subfolder):
    return load_helper(subfolder, 'true_counts.csv',0)

def load_classification(subfolder):
    return load_helper(subfolder, 'group_data.csv',0)

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
    
# run kmeans n times and compute ARI
def kmeans_ARI(recon_matrix, labels_true, num_K, n):
    dist = np.zeros(n)
    kmeans = KMeans(n_clusters=num_K)
    for i in range(n):
        labels_pred = kmeans.fit_predict(recon_matrix)
        dist[i] = adjusted_rand_score(labels_true, labels_pred)
    [lower_quartile, median, upper_quartile] = percentile(dist, [25, 50, 75])
    # calculate min/max
    data_min, data_max = dist.min(), dist.max()
    
    return data_min, lower_quartile, median, upper_quartile, data_max
    
    
    
"""
Here row is cells column are genes
"""
def main():
    report_name = "ARI_report_csv"
    subfolders = ['2_groups_1000_cells_5000_genes/',
                  '5_groups_1000_cells_5000_genes/',
                  '10_groups_1000_cells_5000_genes/',
                  '2_groups_10000_cells_1000_genes/',
                  '5_groups_10000_cells_1000_genes/',
                  '10_groups_10000_cells_1000_genes/'];
#    file_names = ['counts.csv', 'true_counts.csv','output_ALRA_hack.csv','output_sI_thresh.csv']
    file_names = ['output_ALRA_hack.csv','output_sI_thresh.csv']

    num_tries = 50 # how many times we run kmeans
    pca_dim = 2
    pca = PCA(n_components = pca_dim)
    
    data_path = os.path.dirname(os.getcwd())
    data_path = data_path + '/SplatGenData/'
    for subfolder in subfolders:
        true_labels = np.array(load_classification(subfolder))
        true_labels = true_labels.flatten()
        print(np.shape(true_labels))
        num_K = np.max(true_labels) + 1
        results = np.zeros([len(file_names), 5])
        data = {'ReconMatrix': file_names}
        header_lines = 0
        for i in range(len(file_names)):
            print(i)
            if i > 1:
                header_lines = 1
                print(header_lines)
            recon_matrix = load_helper(subfolder, file_names[i],header_lines)
            X_pc = pca.fit_transform(recon_matrix.T)
            results[i,:] = kmeans_ARI(X_pc, 
                                      true_labels, num_K, num_tries)
        
        data["min"] = results[:, 0]
        data["25%"] = results[:, 1]
        data["median"] = results[:,2]
        data["75%"] = results[:,3]
        data["max"] = results[:,4]
        
        df = pd.DataFrame(data)
        df.to_csv(data_path + subfolder + report_name)
        print("finished writing report for", subfolder)
            
        
        
    

if __name__ == '__main__':
    main()
