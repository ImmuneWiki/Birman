import os
import numpy as np
from bam_tools import *
from report_tools import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("bamfile", help="BAM file after comparison")
parser.add_argument("referencefile", help="IGHV on IMGT website")
parser.add_argument("-c", "--CDR", help="draw only CDR area", action="store_true")
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

if not os.path.exists('MutCounts/'):
    os.makedirs('MutCounts/')  # MutCounts文件夹存放突变碱基信息


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

    for fastqseq in fastqseq_list:  # fastqseq_list存放着处理过的fastq序列
        for i in range(0, len(fastqseq)):  # 更新nt_list
            if fastqseq[i] == 'A':
                nt_list[i, 0] += 1
            elif fastqseq[i] == 'T':
                nt_list[i, 1] += 1
            elif fastqseq[i] == 'C':
                nt_list[i, 2] += 1
            elif fastqseq[i] == 'G':
                nt_list[i, 3] += 1
    nt_list = nt_list.astype(int)
    gene_name = gene_header.split('|')[1].replace('*', '_').replace('/', '.')
    if sum(sum(nt_list)) == 0:
        pass
    else:
        CDR_region = find_CDR()
        freqs, depth = getMutFreq(nt_list, reference_gene)
        Mutcounts = freqs * np.transpose([depth, depth, depth, depth])
        Mutcounts = Mutcounts.astype(int)
        np.savetxt('MutCounts/' + gene_name + '.txt', Mutcounts, fmt='%s')
        if args.CDR:
            plotMutations(freqs, gene_name, CDR_region[gene_name])
        else:
            plotMutations(freqs, gene_name)


referencefile.close()
bamfile.close()