import os
import matplotlib.pyplot as plt
import numpy as np
from bam_tools import *
from time import *
import pandas as pd
import re


def plotMutations(Freqs, figName, region=0):
    ### This function returns NumPy arrays of prediction intervals.
    ### Parameters:    Test: Numpy array [length(gene) x 4]
    ###                        Alternative nucleotide frequencies in sample as returned by getAltFreq
    ###                Testd: Numpy array [length(gene) x 1]
    ###                        Position-wise depth in sample as returned by getAltFreq
    ###                gene: string
    ###                        File "gene.fasta" must be in nucleotideCounts directory
    ###                region: range
    ###                        Specify the region of gene sequence.
    ###                figName: string
    ###                        File name to save generated figures

    # Mutcounts = Freqs * np.transpose([Depths, Depths, Depths, Depths])
    colors = ['#57bddb', '#ea4848', '#3bb971', '#f39c12']
    Freqs[0:region[0],:]=[0, 0, 0, 0]
    Freqs[region[1]:region[2],:]=[0, 0, 0, 0]
    Freqs[region[3]:len(Freqs), :] = [0, 0, 0, 0]
    # Plot position-wise mismatch frequency
    plt.plot(Freqs[:, 0] + Freqs[:, 1] + Freqs[:, 2] + Freqs[:, 3], 'b')
    plt.rcParams.update({'font.size': 10})
    plt.axis([0, len(Freqs), 0, 1])
    plt.xlabel('Position (bp)')
    plt.ylabel('Mutation frequency')
    plt.savefig('MutCounts/'+figName + '.pdf', bbox_inches='tight')
    plt.close()

    # Plot position-wise mismatch frequency for each nucleotide
    ymax = 1  # 突变频率上限设定为1
    f, axarr = plt.subplots(2, 2)
    plt.rcParams.update({'font.size': 10})
    axarr[0, 0].plot(Freqs[:, 0], color=colors[0])  # 第0列是A
    axarr[0, 0].set_ylim([0, ymax])
    axarr[0, 0].set_xlim([0, len(Freqs)])
    axarr[0, 0].set_title('Mut. to A',color=colors[0])
    axarr[0, 1].plot(Freqs[:, 1], color=colors[1])  # 第1列是T
    axarr[0, 1].set_ylim([0, ymax])
    axarr[0, 1].set_xlim([0, len(Freqs)])
    axarr[0, 1].set_title('Mut. to T',color=colors[1])
    axarr[1, 0].plot(Freqs[:, 2], color=colors[2])  # 第2列是C
    axarr[1, 0].set_ylim([0, ymax])
    axarr[1, 0].set_xlim([0, len(Freqs)])
    axarr[1, 0].set_title('Mut. to C',color=colors[2])
    axarr[1, 1].plot(Freqs[:, 3], color=colors[3])  # 第3列是G
    axarr[1, 1].set_ylim([0, ymax])
    axarr[1, 1].set_xlim([0, len(Freqs)])
    axarr[1, 1].set_title('Mut. to G',color=colors[3])
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    axarr[1, 1].set_xlabel('Position (bp)')
    axarr[1, 0].set_xlabel('Position (bp)')
    axarr[1, 0].set_ylabel('Mut. frequency')
    axarr[0, 0].set_ylabel('Mut. frequency')
    plt.savefig('MutCounts/'+figName + '_bases.pdf', bbox_inches='tight')
    plt.close()
    return


def distance(seq1, seq2):
    """
    计算距离，忽视了N与ATCG的区别，仅计算了ATCG之间不同的个数，如果情况有变可以进一步修改
    :param seq1: 核苷酸序列1
    :param seq2: 核苷酸序列2
    :return: 碱基突变的坐标（格式为数字+A/T/C/G，表明几号位突变成了什么）
    """
    different_index = []
    different = []
    seq1_list = list(seq1)
    seq2_list = list(seq2)
    for i in range(len(seq2)):
        if seq1_list[i] == 'N' or seq2_list[i] == 'N':
            continue
        if seq1_list[i] != seq2_list[i]:
            different_index.append(i)
            different.append(str(i) + seq2_list[i])
    return different_index, different


def getDistanceMatrix(mutant_index_series, mutant_series, mutant_len_series):
    """
    得到距离矩阵，为了降低时间消耗，所以只计算差别个数之差在1之内的序列，因为差别个数在2以上的序列至少有两个位置不一样，距离至少为2
    :param mutant_index_series:一个Series，包含unique_seq处理过后的['different_index']列
    :param mutant_series:一个Series，包含unique_seq处理过后的['different']列
    :param mutant_len_series:一个Series，包含unique_seq处理过后的['different_indexlen']列
    :return:距离矩阵
    """
    index_dict = {}  # 这个主要是为了降低时间消耗，算出差别在1以内的下一个序列最大的坐标
    mutant_num = list(mutant_len_series)
    for i in mutant_num:
        index_dict[i] = index_dict.get(i, 0) + 1
    key_list = list(index_dict.keys())
    for i in range(len(key_list)):
        if i == 0:
            continue
        index_dict[key_list[i]] = index_dict[key_list[i - 1]] + index_dict[key_list[i]]
    for i in range(len(key_list) - 1):
        index_dict[key_list[i]] = index_dict[key_list[i + 1]]


    LenRow = mutant_index_series.shape[0]
    Dis_Mat = np.zeros((LenRow, LenRow))
    Mut_Mat = np.zeros((LenRow, LenRow))
    for i in range(0, LenRow):
        for j in range(i+1, index_dict[mutant_num[i]]):  # 为了计算出距离，需要同时用到index信息与mutant信息
            all_index = len(set(mutant_index_series[i] + mutant_index_series[j]))  # 去重的index数量
            same_mutant = len(mutant_series[i] + mutant_series[j]) - len(set(mutant_index_series[i] + mutant_index_series[j]))  # 相同的mutant数量
            diff = all_index - same_mutant  # 两者相减即为距离
            Dis_Mat[i, j] = diff
            Mut_Mat[i, j] = same_mutant  # 两个列表都存在的元素，即共同突变
    return Dis_Mat, Mut_Mat



if __name__ == "__main__":
    begin_time = time()
    merge_seq_df = pd.read_csv("merge_seq.csv")
    seqs = merge_seq_df['merge_align_seq']
    merge_seq_df.drop(['merge_align_quality'], axis=1, inplace=True)
    reference_seq = 'CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCC'
    different_list = []
    different_index_list = []
    different_listlen = []
    different_str_list = []
    for seq in seqs:
        different_index, different = distance(reference_seq, seq)  # different_index是一个列表，记录着这个序列的突变位点的坐标和类型
        different_str_list.append(''.join(different))
        different_list.append(different)
        different_index_list.append(different_index)
        different_listlen.append(len(different))
    merge_seq_df['different_str'] = different_str_list
    merge_seq_df['different'] = different_list
    merge_seq_df['different_index'] = different_index_list
    merge_seq_df['different_indexlen'] = different_listlen

    unique_seq_df = merge_seq_df.merge(merge_seq_df.groupby(['different_str']).size().reset_index().rename(columns={0: 'frequency'}))
    unique_seq_df.loc[unique_seq_df['different_str'].duplicated(), 'frequency'] = -1
    unique_seq_df.drop(unique_seq_df[unique_seq_df.frequency == -1].index, inplace=True)
    unique_seq_df.reset_index(drop=True, inplace=True)
    # unique_seq_df.to_csv("unique_seq.csv")
    mid_time = time() - begin_time
    print('该循环程序运行时间：%.8f' % mid_time)
    print(unique_seq_df)
    # unique_seq_df = pd.read_csv('unique_seq.csv')
    # 得到BCR间距离矩阵
    unique_seq_df.sort_values(by='different_indexlen', ascending=True, inplace=True, ignore_index=True)  # 按照突变个数从小到大排序
    print(unique_seq_df)
    # unique_seq_df.to_csv("unique_seq.csv")
    DistanceMatrix, Mut_same_Matrix = getDistanceMatrix(unique_seq_df['different_index'], unique_seq_df['different'], unique_seq_df['different_indexlen'])  # 计算距离矩阵与突变相似性矩阵
    DistanceMatrix = np.mat(DistanceMatrix)
    np.savetxt('DistanceMatrix.txt', DistanceMatrix, fmt='%d')



    # 得到突变相似性矩阵
    Mutant_count = np.array(unique_seq_df['different_indexlen'])
    Mut_same_Matrix = np.mat(Mut_same_Matrix)
    Mutant_count = np.maximum(Mutant_count, 1)  # 把突变数量为0的那个序列置1，因为下面一步要做分母
    Mut_prop_Matrix = (Mut_same_Matrix + Mut_same_Matrix.T) / Mutant_count
    np.savetxt('Mutant_prop_Matrix.txt', Mut_prop_Matrix, fmt='%.3f')


