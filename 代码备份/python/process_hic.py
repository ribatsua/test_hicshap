#由hic文件提取稀疏矩阵按染色体保存整个数据集
import hicstraw
import os
import numpy as np
from scipy.sparse import csr_matrix
from tqdm import tqdm
import math

#计算每条染色体在分辨率res下划分的bin的数量
def create_bed(res,genome_reference_path):
    file_path=genome_reference_path
    # file_path = '/home/python/higashi/cellcycle/250k/config/mm9.chrom.sizes.txt'  
    with open(file_path, 'r') as file:
        lines = file.readlines()
    chromosome_bins = {}
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
    return chromosome_bins
#按文件名对文件列表排序
def sort_key(path):
    # 先使用斜杠分割路径，然后使用下划线分割每个部分，获取数字作为排序关键字
    parts = []
    for part in path.split('/'):
        parts.extend(part.split('_'))
    return [int(part) if part.isdigit() else part for part in parts]
#获取数据集根目录下所以chrom命名的txt文件并返回排序后的文件列表
def find_file(folder_path):
    file_path=[]
    for root, dirs, files in os.walk(folder_path):
            for file in files:
                # if file.endswith('.txt') and chrom == file.split(".")[0].split("_")[-1]: #添加.split("_")[-1]后通用于以cell_1_chr1.txt命名文件
                if file.endswith('.hic') and 'cortex' == file.split(".")[0].split("-")[0].split("_")[-1]:
                # if file.endswith('.txt'):
                    file_path.append(os.path.join(root, file))
    sorted_paths = sorted(file_path, key=sort_key)
    return sorted_paths          

def process_chrom(chrom,res,sorted_paths,bin_count,output_dir):
    sparse_matrices_list=[]
    chrNum = chrom.split("chr")[1]
    for path in tqdm(sorted_paths,desc="Processing %s"%chrom):
        result = hicstraw.straw("observed", 'NONE', path, chrNum,chrNum, 'BP', res)
        sorted_result = sorted(result, key=lambda x: (x.binX, x.binY))
        rows = []
        cols = []
        values = []
        for elem in sorted_result:
            row, col, value = (int(elem.binX/res),int(elem.binY/res),elem.counts)
            rows.append(row)
            cols.append(col)
            values.append(value)
        sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(bin_count, bin_count), dtype=np.float32)#统一稀疏矩阵形状为bin数量
        sparse_matrices_list.append(sparse_matrix_csr)

    np.save(os.path.join(output_dir,'%s_sparse_matrices.npy'% chrom), sparse_matrices_list)
    print (chrom, "finished")

def main():
    genome_reference_path='/home/python/higashi/dataset2/cortex/config/mm10.chrom.sizes.txt'
    folder_path = '/home/python/higashi/dataset_hic/dataset2/hic_data'
    output_dir = '/home/python/higashi/dataset_hic/dataset2/cortex50k/raw'

    res = 50000
    chrom_list = ["chr1","chr2","chr3","chr4","chr5",
    "chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15",
    "chr16","chr17","chr18","chr19","chrX"]
    chrom_list = np.array(chrom_list)
    print(chrom_list)
    sorted_paths=find_file(folder_path)
    chromosome_bins = create_bed(res,genome_reference_path)
    for chrom in chrom_list:
        bin_count=chromosome_bins[chrom]
        process_chrom(chrom,res,sorted_paths,bin_count,output_dir)

if __name__ == '__main__':
    main()


    





        




##############################################dchic######################################################

# Extract Data
# file_path='/home/python/higashi/dataset_hic/dataset2/hic_data/GSM4382149_cortex-p001-cb_001.contacts.hic'
# hic = hicstraw.HiCFile(file_path)
# genome = hic.getGenomeID()
# resolutions = hic.getResolutions()

# print(f" - These are the resolutions of your file: {resolutions}")
# print(f" - This is the genome of your file: {genome}.")

# # for chrom in hic.getChromosomes():
# #     print(chrom.name)
# chrs = []
# sizes = []
# chrSizes = {}

# # # Removing Chromosomes 
# remove='Y,M'
# chrsToRemove = []
# if remove:
#     chrsToRemove = remove.strip().lower().split(",")
# chrsToRemove.append("all") # default
# chrsToRemove.append("mt")
# print('chrsToRemove:',chrsToRemove)
# print(" - Removing these chromosomes:")
# for chrom in hic.getChromosomes():
#     if chrom.name.lower() in chrsToRemove:
#         if chrom.name.lower() == "all":
#             print(f"  - {chrom.name} (.hic file artifact removed by default)")
#         else:
#             print(f"  - {chrom.name}")
#     else:
#         chrTag = "chr" + chrom.name if "chr" not in chrom.name else chrom.name
#         chrSize = chrom.length
#         chrSizes[chrTag] = chrSize
#         chrs.append(chrTag)
#         sizes.append(chrSize)
# print(chrs)
# # Make bed files

# # name = results.prefix + "_" + str(results.res) + "_abs.bed" #data_200000_abs.bed
# iterator = 1 # one-based
# res = 250000

# positionHash = {} # list of hash tables (first coord : index), one per chromosome

# # with open(name, "w") as bedfile:
##############################################dchic#####################################################
# print(" - Creating bed file")
# for a in range(len(chrs)):
#     chrDic = {}
#     length = chrSizes.get(chrs[a])
#     posIterator = int(res)
#     while posIterator <= length:
#         chrDic[(posIterator-res)] = iterator
#         # bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(posIterator) + "\t" + str(iterator) + "\n")
#         posIterator+=res
#         iterator+=1
#     if (posIterator-res) < length:
#         chrDic[(posIterator-res)] = iterator
#         # bedfile.write(chrs[a] + "\t" + (str(posIterator-res)) + "\t" + str(length) + "\t" + str(iterator) + "\n")
#         iterator +=1
#     positionHash[chrs[a]] = chrDic
# print(positionHash[chrs[0]])
##############################################dchic#####################################################

####################################cp##################################################################
###计算每个染色体划分的bin数量，用于后续指定创建稀疏矩阵形状
# import math
# chromosome_bins = {}
# for a in chrs:
#     length = chrSizes.get(a)
#     bin_count = math.ceil(int(length) / res)
#     chromosome_bins[a]=bin_count
# # print(chromosome_bins['chr1'])
####################################cp##################################################################


####################################cp##################################################################
#     print(" - Bed File Creation Done.")

# # Sparse Matrix Dump

# outname = results.prefix + "_" + results.res + ".matrix"
# outname = 'test01' + "_" + '250k' + ".matrix"

# chrList = []
# print(" - Processing These Chromosomes: ")
# for c in range(len(chrs)):
#     chr = chrs[c]
#     chrSpecific = []
#     chrNum = chr.split("chr")[1]
#     print(f"  - {chr}")
#     result = hicstraw.straw("observed", 'NONE', results.file, chrNum, chrNum, 'BP', res)
#     for elem in result:
#         chrSpecific.append(elem)
#     chrSpecific.sort(key = lambda x : (x.binX, x.binY))
#     chrList.append(chrSpecific)

# with open(outname, "w") as outfile:
#     for c in range(len(chrList)):
#         chr = chrs[c]
#         for elem in chrList[c]:
#             outfile.write("{0}\t{1}\t{2}\n".format(positionHash[chr][elem.binX], positionHash[chr][elem.binY], int(elem.counts)))


######################################################################################################################

# chr='chr1'
# print(f"  - {chr}")
# result = hicstraw.straw("observed", 'NONE', file_path, '1','1', 'BP', 250000)
# chr_path='/home/python/higashi/dataset_hic/python/test.chr1.txt'
# # for i in range(len(result)):
#     # print("{0}\t{1}\t{2}".format(result[i].binX/250000, result[i].binY/250000, result[i].counts))
# # 按照 binX 和 binY 属性进行排序
# sorted_result = sorted(result, key=lambda x: (x.binX, x.binY))

######################################################################################################################


# # 将排序后的结果写入文件
# with open(chr_path, 'w') as file:
#     for elem in sorted_result:
#         file.write(f'{int(elem.binX/250000)}\t{int(elem.binY/250000)}\t{elem.counts}\n')
##########################################cp_test######################################################
# import numpy as np
# from scipy.sparse import csr_matrix

# path='/home/python/higashi/dataset_hic/python/test.chr1.txt'
# sparse_matrices_list=[]
# bin_count=chromosome_bins['chr1']
# with open(path, 'r') as f:
#     lines = f.readlines()
#     rows = []
#     cols = []
#     values = []
#     for line in lines:
#         # row, col, value = map(float, line.split())
#         row, col, value = map(lambda x: int(x) if '.' not in x else float(x), line.split('\t'))#将行列索引转换成整数
#         rows.append(row)
#         cols.append(col)
#         values.append(value)
#     # sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(max(rows)+1, max(cols)+1), dtype=np.float32)
#     sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(bin_count, bin_count), dtype=np.float32)#统一稀疏矩阵形状为bin数量
#     sparse_matrices_list.append(sparse_matrix_csr)
# print(sparse_matrices_list)
########################################cp_test################################################################

##########################################cp_test######################################################
# import numpy as np
# from scipy.sparse import csr_matrix

# path='/home/python/higashi/dataset_hic/python/test.chr1.txt'

#测试从.hic中提取的数据直接转换成稀疏矩阵
# chr='chr1'
# sparse_matrices_list=[]
# bin_count=chromosome_bins[chr]
# chrNum = chr.split("chr")[1]
# print(f"  - {chr}")
# result = hicstraw.straw("observed", 'NONE', file_path, chrNum,chrNum, 'BP', 250000)
# # chr_path='/home/python/higashi/dataset_hic/python/test.chr1.txt'
# # for i in range(len(result)):
#     # print("{0}\t{1}\t{2}".format(result[i].binX/250000, result[i].binY/250000, result[i].counts))
# # 按照 binX 和 binY 属性进行排序
# sorted_result = sorted(result, key=lambda x: (x.binX, x.binY))

#     # file.write(f'{int(elem.binX/250000)}\t{int(elem.binY/250000)}\t{elem.counts}\n')
#     # lines = f.readlines()
# rows = []
# cols = []
# values = []
# for elem in sorted_result:
#     # row, col, value = map(float, line.split())
#     row, col, value = (int(elem.binX/250000),int(elem.binY/250000),elem.counts)
#     rows.append(row)
#     cols.append(col)
#     values.append(value)
# # sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(max(rows)+1, max(cols)+1), dtype=np.float32)
# sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(bin_count, bin_count), dtype=np.float32)#统一稀疏矩阵形状为bin数量
# sparse_matrices_list.append(sparse_matrix_csr)
# print(sparse_matrices_list[0])
########################################cp_test################################################################



# with open(chr_path, 'w') as file:
#     for i in range(len(result)):
#             file.write(f'{int(result[i].binX/250000)}\t{int(result[i].binY/250000)}\t{result[i].counts}\n')
