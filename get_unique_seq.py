from report_tools import *
from bam_tools import *
import os
import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("bamfile", help="BAM file after comparison")
parser.add_argument("referencefile", help="IGHV on IMGT website")
parser.add_argument("-b", "--txt", help="the BAM file was viewed", action="store_true")
args = parser.parse_args()
if args.txt:
    bamfile = open(args.bamfile)
else:
    os.system("samtools view "+ args.bamfile + " > " + args.bamfile + ".txt")
    bamfile = open(args.bamfile.txt)
referencefile = open(args.referencefile)


reference_genelist = referencefile.readlines()
line = bamfile.readline()
line_list = line.split('\t')

if not os.path.exists('unique_seq/'):
    os.makedirs('unique_seq/')  # unique_seq文件夹
if not os.path.exists('unique_seq_pkl/'):
    os.makedirs('unique_seq_pkl/') # unique_seq_pkl文件夹

for i in range(0, len(reference_genelist) // 2):  # 2i是基因名，2i+1是基因
    fastqseq_list = []
    quality_list = []
    seq_id_list = []
    gene_header = reference_genelist[2 * i].strip()
    reference_gene = reference_genelist[2 * i + 1].strip()
    nt_list = np.zeros((len(reference_gene), 4))  # nt_list用于存放碱基统计信息，参考基因组的长度行，4列，行索引代表碱基位置，列索引代表ATCG
    while True:
        if not line:  # 读到末尾停止
            break
        if line_list[2] not in gene_header:  # 注意一定要输入“排序”后的bam文件，因为这样的文件是把同样的BCR放一起了，然后这句话的意思就是读到这个BCR读段末尾就停止
            break
        if line_list[4] == 0:  # 去除匹配质量为0的序列
            line = bamfile.readline()
            line_list = line.split('\t')
            continue

        fastqseq, quality, seq_id = align(reference_gene, line_list[9], line_list[3], line_list[5], line_list[10], line_list[0])
        fastqseq_list.append(fastqseq)
        quality_list.append(quality)
        seq_id_list.append(seq_id)

        line = bamfile.readline()
        line_list = line.split('\t')

    gene_name = gene_header.split('|')[1].replace('*', '_').replace('/', '.')

    seq_df = pd.DataFrame({'seq_id': seq_id_list,
                       'align_seq': fastqseq_list,
                       'align_quality': quality_list})
    merge_seq_df = pairs_merge(seq_df)
    seqs = merge_seq_df['merge_align_seq']
    merge_seq_df.drop(['merge_align_quality'], axis=1, inplace=True)
    different_list = []
    different_index_list = []
    different_listlen = []
    different_str_list = []
    for seq in seqs:
        different_index, different = distance(reference_gene, seq)  # different_index是一个列表，记录着这个序列的突变位点的坐标和类型
        different_str_list.append(''.join(different))
        different_list.append(different)
        different_index_list.append(different_index)
        different_listlen.append(len(different))
    merge_seq_df['different_str'] = different_str_list
    merge_seq_df['different'] = different_list
    merge_seq_df['different_index'] = different_index_list
    merge_seq_df['different_indexlen'] = different_listlen

    unique_seq_df = merge_seq_df.merge(
        merge_seq_df.groupby(['different_str']).size().reset_index().rename(columns={0: 'frequency'}))
    unique_seq_df.loc[unique_seq_df['different_str'].duplicated(), 'frequency'] = -1
    unique_seq_df.drop(unique_seq_df[unique_seq_df.frequency == -1].index, inplace=True)
    unique_seq_df.reset_index(drop=True, inplace=True)
    unique_seq_df.to_pickle('unique_seq_pkl/' + gene_name + 'unique_seq.pkl')
    unique_seq_df.to_csv('unique_seq/' + gene_name + 'unique_seq.csv')

referencefile.close()
bamfile.close()

