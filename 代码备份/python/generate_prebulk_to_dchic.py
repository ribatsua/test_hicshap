############################################################
####
#####由原始的稀疏矩阵生成prebulk数据及bed文件作为dchic的输入
#####
##########################################################
#############################################################
#############################################################
##################################################################
##################################################################
import os
import numpy as np
import pandas as pd
import scipy.sparse
import math

from tqdm import tqdm


################从基因组参考文件计算一定分辨率下各染色体bin的数量
def create_bed(res,genome_reference_path):
    file_path=genome_reference_path 
    with open(file_path, 'r') as file:
        lines = file.readlines()
    chromosome_bins = {}
    max_lens = {}
    # 特定分辨率
    resolution = res # 请替换成你的分辨率
    # 遍历文件中的每一行
    for line in lines:
        # 分割每一行的染色体和对应的长度
        chromosome, length = line.strip().split('\t')
        
        # 计算染色体对应的bin数量
        # bin_count = int(length) // resolution
        bin_count = math.ceil(int(length) / resolution)
        
        # 存储结果到字典中
        chromosome_bins[chromosome] = bin_count
        max_lens[chromosome] = int(length)
    # bin_chrom_list = []
    # bin_start_list = []
    # bin_end_list = []
    # for chrom in chromosome_bins:
    #     bin_count = chromosome_bins[chrom]
    #     bin_chrom_list+=[chrom]*bin_count
    #     bin_start_list.extend((np.arange(bin_count) * res))
    #     end=(np.arange(bin_count) + 1) * res
    #     # print(len(end))
    #     max_len=max_lens[chrom]
    #     end[bin_count-1]=max_len
    #     # bin_end_list.extend(((np.arange(bin_count) + 1) * res))
    #     bin_end_list.extend(end)
    #     # print (chrom, "finished")
    #     # print(chrom)
    #     # print(bin_count)
    # index_list=np.arange(1,len(bin_chrom_list)+1)
    # df=pd.DataFrame()
    # df['chr']=bin_chrom_list
    # df['start']=bin_start_list
    # df['end']=bin_end_list
    # df['index']=index_list
    # # bed_path='test_bed.bed'
    # bed_path=os.path.join(bed_dir, 'test.bed')
    # df.to_csv(bed_path, sep='\t', index=None, header=False, mode='w')
    # print(f"bed file written to {bed_path}")
    return chromosome_bins


def process_chrom(chrom,cell_type):
# Get the raw sparse mtx list
    origin_sparse = np.load(os.path.join(temp_dir, "%s_sparse_matrices.npy" % chrom), allow_pickle=True)
    size = origin_sparse[0].shape[0]
    pre_bulk_list=[]
    # sum_by_label=[]
    a={}
    b={}
    # c={}
    for j in list(set(cell_type)):
        # print(j)
        indices = [index for index, value in enumerate(cell_type) if value == j]
        a[j]=indices
        b[j] = np.zeros_like(np.array(origin_sparse[0]))
        for i in tqdm(a[j]):
            proba = np.array(origin_sparse[i])    
            b[j] +=proba
        temp = np.array(b[j].item().todense())
        temp=scipy.sparse.csr_matrix(temp)
            # temp = b[j]
        pre_bulk_list.append(temp)
    return pre_bulk_list



import pickle
res=50000
# temp_dir='/home/python/higashi/cellcycle/250k/raw'
temp_dir = '/home/python/higashi/dataset_hic/dataset2/cortex50k/raw'
output_dir='/home/python/higashi/dcpc/from_hicfile_cortex50k_mat'
# bed_dir='/home/python/higashi/dcpc/cellcycle250k/Data'
genome_reference_path = '/home/python/higashi/dataset2/cortex/config/mm10.chrom.sizes.txt'
if not os.path.exists(output_dir):
    # 如果路径不存在，创建路径
    os.makedirs(output_dir)
    print(f"目录 '{output_dir}' 已创建.")
else:
    print(f"目录 '{output_dir}' 已存在.")
chrom_list = ["chr1","chr2","chr3","chr4","chr5",
    "chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15",
    "chr16","chr17","chr18","chr19"]

# chrom_list = ['chr1', 'chr2']

chrom_list = np.array(chrom_list)
print('chrom_list:',chrom_list)
# with open('matrix_values.txt', 'w') as file:
# cell_type_info=pickle.load(open(os.path.join(data_dir, "label_info.pickle"), "rb"))
#使用重复样本标签

# cell_type_info=pickle.load(open('/home/python/higashi/cellcycle/250k/config/label_info.pickle', "rb"))
cell_type_info=pickle.load(open('/home/python/higashi/dataset_hic/dataset2/cortex250k/label_info.pickle', "rb"))

cell_type=cell_type_info['cell_type']
cell_type_list=[]
chromosome_bins=create_bed(res,genome_reference_path)
for j in list(set(cell_type)):
    cell_type_list.append(j)
bin_count=0
pre_bulk_raw=[]
for chrom in chrom_list:
    pre_bulk_list=process_chrom(chrom,cell_type)
    print (chrom, "finished")
    pre_bulk_list_unzip=[]
    for bulk in pre_bulk_list:
        rows, cols, values = scipy.sparse.find(bulk)
        rows+=bin_count
        cols+=bin_count
        pre_bulk_list_unzip.append([cols,rows,values])
    bin_count += chromosome_bins[chrom] 
    pre_bulk_raw.append(pre_bulk_list_unzip)
temp_bulk ={}
for i in tqdm(range(len(pre_bulk_raw[0]))):
    sample_name=cell_type_list[i]+'.matrix'
    save_path=os.path.join(output_dir, sample_name)
    # save_path=cell_type_list[i]+'_prebulk_mat.txt'   
    temp_bulk[i]=[]
    for j in range(len(pre_bulk_raw)):
        # print('j,i',j,i)
        temp_bulk[i].append(pre_bulk_raw[j][i])
    temp_bulk[i]=np.column_stack([np.vstack(group) for group in temp_bulk[i]])
    np.savetxt(save_path, temp_bulk[i].T, delimiter='\t', fmt='%d')
    print(f"Data written to {save_path}")