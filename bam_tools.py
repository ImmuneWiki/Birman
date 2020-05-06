import re
import numpy as np
import pandas as pd
from report_tools import *
import time
import sys

def align(reference_seq, seq, pos, CIGAR, quality, seq_id):
    """
    以bam文件提供的参数进行联配，为了显示CDR区，将未测序部分的碱基填充N，碱基质量填充！（即最低质量）
    :param reference_seq: 并未用到，只不过留着说不定会有用，没用可以直接忽视
    :param seq: bam文件的第十列，核苷酸序列
    :param pos: bam文件的第四列，开始比对上的位置
    :param CIGAR: bam文件的第六列，比对信息，S代表去除的碱基个数，M代表匹配/错配的碱基个数，D代表缺失的碱基个数，I代表插入的碱基个数（因为插入不好处理，我忽视了插入的情况，可以以后优化）
    :param quality: bam文件的第十一列，碱基质量
    :param seq_id: bam问价的第一列，读段id
    :return: (填充N后的序列,填充!后的碱基质量，读段id)
    """
    CIGAR_list = re.findall(r'[0-9]+|[A-Z]+', CIGAR)
    pos = int(pos)
    seq_new = ''
    quality_new = ''
    for i in range(0, len(CIGAR_list) // 2):
        if CIGAR_list[i*2+1] == 'S':
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[2*i+1] == 'M':
            seq_new = seq_new + seq[0:int(CIGAR_list[2*i])]
            quality_new = quality_new + quality[0:int(CIGAR_list[2*i])]
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[2*i+1] == 'D':
            seq_new = seq_new + int(CIGAR_list[2*i]) * 'N'
            quality_new = quality_new + int(CIGAR_list[2*i]) * '!'
        elif CIGAR_list[2*i+1] == 'I':
            seq = seq[int(CIGAR_list[2*i]):]
            quality = quality[int(CIGAR_list[2*i]):]
        elif CIGAR_list[len(CIGAR_list)-1] == 'S':
            pass


    seq_new = (pos-1)*'N'+seq_new  # 在序列前填充N
    seq_new = seq_new.ljust(len(reference_seq), 'N')  # 在序列后填充N
    quality_new = (pos-1)*'!'+quality_new  # 在序列前填充!
    quality_new = quality_new.ljust(len(reference_seq), '!')  # 在序列后填充!
    return (seq_new, quality_new, seq_id)



def getMutFreq(Counts, reference_gene):
    """
    统计各个位点的突变的碱基频率，输入的是nt_list，四列n行的数组
    :param sample: 样品名称
    :return: 频率（n行4列的数组）与深度（一位数组）
    """
    # Compute nucleotide frequencies and depth
    jianji = "ATCGN"
    char_to_int = dict((c, i) for i, c in enumerate(jianji))
    reference_gene_Count = np.array([char_to_int[nt] for nt in reference_gene])
    reference_gene_Count = np.eye(5)[reference_gene_Count.reshape(-1)]  # 对参考序列进行独热编码
    Depths = []
    Freqs = []
    for count in Counts:
        Depths.append(sum(count))
        Freqs.append([c / sum(count) if sum(count) > 0 else 1 for c in count])
    Freqs = np.array(Freqs)
    (n, m) = np.shape(Freqs)
    for i in range(n):
        for j in range(m):
            if reference_gene_Count[i, j] == 1:
                Freqs[i, j] = 0
    return (Freqs, Depths)


def find_CDR(filename="Human IGHV F+ORF+in-frame P含..txt"):
    """
    得到不同BCR的CDR区位置
    :param filename: IMGT下载的填充.的原始BCR序列
    :return: 一个字典，键为BCR序列名称，值为CDR1、CDR2的起止位置
    """
    referencefile = open(filename)
    reference_genelist = referencefile.readlines()
    CDR_dict = {}
    for i in range(0, len(reference_genelist) // 2):  # 2i是基因名，2i+1是基因
        CDR1_start = 27
        CDR1_end = 113
        CDR2_start = 167
        CDR2_end = 194
        gene_header = reference_genelist[2 * i].strip()
        reference_gene = reference_genelist[2 * i + 1].strip()
        gene_name = gene_header.split('|')[1].replace('*', '_').replace('/', '.')
        for i in range(len(reference_gene)):
            if reference_gene[i] == '.' and i < CDR1_start:
                CDR1_start = CDR1_start - 1
                CDR1_end = CDR1_end - 1
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i >= CDR1_start and i <= CDR1_end:
                CDR1_end = CDR1_end - 1
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i < CDR2_start:
                CDR2_start = CDR2_start - 1
                CDR2_end = CDR2_end - 1
            elif reference_gene[i] == '.' and i >= CDR2_start and i <= CDR2_end:
                CDR2_end = CDR2_end - 1
        CDR_dict[gene_name] = [CDR1_start,CDR1_end,CDR2_start,CDR2_end]
    return CDR_dict


def pairs_merge(row_seqdatafame):
    """
    双端测序的结果含有序列相同的读段，本函数是将序列相同的读段合并，如果重合区域碱基不同，则选质量高的碱基
    :param row_seqdatafame: 经过align操作过后的读段dataframe，三列分别为：'seq_id'，'merge_align_seq'，'merge_align_quality'
    :return: 一个dataframe，三列分别为：'seq_id'，'merge_align_seq'，'merge_align_quality'
    """
    id_list = list(set(row_seqdatafame['seq_id']))  # 去除重复id
    pair_id_list = []
    seq_merge_list = []
    quality_merge_list = []

    for id in id_list:
        seq_merge = []
        quality_merge = []
        same_seq_index  = row_seqdatafame[row_seqdatafame['seq_id'] == id].index.tolist()
        if len(same_seq_index) == 2:
            index1 = same_seq_index[0]
            index2 = same_seq_index[1]
            seq1 = list(row_seqdatafame['align_seq'][index1])
            seq2 = list(row_seqdatafame['align_seq'][index2])
            quality1 = list(row_seqdatafame['align_quality'][index1])
            quality2 = list(row_seqdatafame['align_quality'][index2])
            length = len(seq1)
            N_count = 0


            for i in range(0, length):
                if quality1[i] or quality2[i] == '!':
                    N_count += 1
                if quality1[i] >= quality2[i]:
                    seq_merge.append(seq1[i])
                    quality_merge.append(quality1[i])
                else:
                    seq_merge.append(seq2[i])
                    quality_merge.append(quality2[i])


            seq_merge = ''.join(seq_merge)
            if not bool(re.search('NNNNNNNNNNNNNNNNNNNN', seq_merge)):
                pair_id_list.append(id)
                seq_merge_list.append(seq_merge)
                quality_merge = ''.join(quality_merge)
                quality_merge_list.append(quality_merge)

    merge_seq_df = pd.DataFrame({'seq_id': pair_id_list,
                                 'merge_align_seq': seq_merge_list,
                                 'merge_align_quality': quality_merge_list})
    return merge_seq_df




if __name__ == "__main__":
    bamfile = open("SRR4242055.1.sort_viewbam.txt")
    referencefile = open("reference.txt")

    reference_genelist = referencefile.readlines()
    line = bamfile.readline()
    line_list = line.split('\t')


    fastqseq_list = []
    quality_list = []
    seq_id_list = []
    gene_header = '>X60503|IGHV1-18*02|Homo sapiens|F|V-REGION|142417|276 nt|1| | | | |276+24=300|partial in 3| |'
    reference_gene = 'CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTAAGATCTGACGACACGGCC'
    nt_list = np.zeros((len(reference_gene), 4))  # nt_list用于存放碱基统计信息，参考基因组的长度行，4列，行索引代表碱基位置，列索引代表ATCG
    while True:
        if line_list[4] == 0:  # 去除匹配质量为0的序列
            line = bamfile.readline()
            line_list = line.split('\t')
            continue
        if line_list[2] in '>M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188483|296 nt|1| | | | |296+24=320| | |':  # 注意一定要输入“排序”后的bam文件，因为这样的文件是把同样的BCR放一起了，然后这句话的意思就是读到这个BCR读段末尾就停止
            line = bamfile.readline()
            line_list = line.split('\t')
            continue
        if line_list[2] in '>X60503|IGHV1-18*02|Homo sapiens|F|V-REGION|142..417|276 nt|1| | | | |276+24=300|partial in 3':
            fastqseq, quality, seq_id = align(reference_gene, line_list[9], line_list[3], line_list[5], line_list[10], line_list[0])
            fastqseq_list.append(fastqseq)
            quality_list.append(quality)
            seq_id_list.append(seq_id)
            line = bamfile.readline()
            line_list = line.split('\t')
        if line_list[2] in '>HM855463|IGHV1-18*03|Homo sapiens|F|V-REGION|21..316|296 nt|1| | | | |296+24=320| | |':
            break

    seq_df = pd.DataFrame({'seq_id': seq_id_list,
                       'align_seq': fastqseq_list,
                       'align_quality': quality_list})
    merge_seq_df = pairs_merge(seq_df)
    merge_seq_df.to_csv('merge_seq.csv', index=False)

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
    unique_seq_df.to_csv('unique_seq.csv')


    referencefile.close()
    bamfile.close()

