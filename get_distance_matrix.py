from bam_tools import *
from report_tools import *
import pandas as pd
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("unique_file_fold", help="the .pkl files folder after the get_unique_seq operation")
args = parser.parse_args()
unique_file_fold = args.unique_file_fold
unique_file_list = os.listdir(unique_file_fold)

if not os.path.exists('DistanceMatrix/'):
    os.makedirs('DistanceMatrix/')

for unique_file in unique_file_list:
    unique_seq_df = pd.read_pickle(unique_file_fold + "/" + unique_file)
    # 得到BCR间距离矩阵
    seqs = unique_seq_df['merge_align_seq']
    unique_seq_df.sort_values(by='different_indexlen', ascending=True, inplace=True, ignore_index=True)  # 按照突变个数从小到大排序
    DistanceMatrix, Mut_same_Matrix = getDistanceMatrix(unique_seq_df['different_index'], unique_seq_df['different'], unique_seq_df['different_indexlen'])  # 计算距离矩阵与突变相似性矩阵
    DistanceMatrix = np.mat(DistanceMatrix)
    np.savetxt('DistanceMatrix/' + unique_file + 'DistanceMatrix.txt', DistanceMatrix, fmt='%d')


    # 得到突变相似性矩阵
    Mutant_count = np.array(unique_seq_df['different_indexlen'])
    Mut_same_Matrix = np.mat(Mut_same_Matrix)
    Mutant_count = np.maximum(Mutant_count, 1)  # 把突变数量为0的那个序列置1，因为下面一步要做分母
    Mut_prop_Matrix = (Mut_same_Matrix + Mut_same_Matrix.T) / Mutant_count
    np.savetxt('DistanceMatrix/' + unique_file + 'Mutant_prop_Matrix.txt', Mut_prop_Matrix, fmt='%.3f')
