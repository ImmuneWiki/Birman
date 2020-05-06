# Birman
Tools for BCR SHM profiling. We provide a tool for analyzing somatic hypermutation and some other BCR analysis tools. These tools are based on the Python interpreter and need to input the bam file obtained by the comparison software(like bowtie2) and samtools(Note that you need to enter the bam file after sort). The program functions are: find SHM, output BCR Sequence and its abundance, drawing mutations, drawing BCR network diagrams, etc.  
我们提供了分析体细胞超突变的一个工具以及其他的一些BCR分析工具，这些工具基于python的解释器，需要输入bowtie2等比对软件处理，以及samtools转化得到的bam文件(注意需要输入sort过后的bam文件)，程序的功能有：找出SHM、输出BCR序列及其丰度、绘制突变情况、绘制BCR网络图等。
## Dependencies
1. Python 3
2. Numpy
3. Pandas
4. Matplotlib
5. Re
6. Networkx(Only for get_networks program)

## List of Tools
- __plot_mutation:__  
  - Description: This tool can identify high mutation sites and plot mutations.  
  - Usage: 
- __get_unique_seq:__  
  - Description: This tool can count sequence types and abundance.
  - Usage: 
- __get_distance_matrix:__  
  - Description: This tool can get the distance (edit distance) between sequences. The distance matrix can be used alone, or it can be used with the fourth program to draw a BCR network diagram to see the mutation.  
此工具可以得到序列间距离（编辑距离），距离矩阵可单独使用，也可配合第四个程序绘制BCR网络图，看突变情况。
  - Usage: 

- __get_networks:__  
  - Description:This tool can get the distance network relationship between sequences. By default, the sequence network with distance 1 is drawn, and the parameter statistical distance network with sequence 2 can be adjusted.  
此工具可以得到序列间距离网络关系，默认绘制距离为1的序列网络，可以调整参数统计距离为2的序列网络。
  - Usage: 





