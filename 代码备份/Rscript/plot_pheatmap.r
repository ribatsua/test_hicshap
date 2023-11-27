library(rhdf5)
library(vegan)
library(pheatmap)
library(data.table)
library(ggplot2) 


method_id_list<-c('manhattan','euclidean','canberra','gower','mahalanobis','chord')

id<-c(seq(4041,4045))
# id<-192
bin_id<-paste0('bin_',id)
bin_id
group_id<-'compartment_raw'
data_dir<-'/home/python/higashi/dataset_hic/dataset2/cortex250k/compartment_raw/'

df<-paste0(data_dir,'ori_shap_values.h5')

file <- H5Fopen(df)
group_info <- H5Gopen(file, group_id)
rhs <- group_info$label$cell_type



#####################################################################################################
# data<-fread("/home/compartment_zscore_shap_hic.csv")
# filtered_data <- subset(data, effect_size>0.16)
# print(paste('filtered bin_num:',length(filtered_data$bin)))
# # length(filtered_data$bin)
# filtered_data <- filtered_data[order(filtered_data$effect_size, decreasing = TRUE), ]
# n <- 100 
# top_n_bins <- head(filtered_data, n)
# bin_id<-top_n_bins$bin


# group_id<-'compartment_zscore'
# df<-paste0('/home/python/higashi/dataset2/cortex/250k/data/',group_id,'/shap_values.h5')
# file <- H5Fopen(df)
# group_info <- H5Gopen(file, group_id)
# rhs <- group_info$label$cell_type

#####################################################################################################


for (j in method_id_list){
        for (i in bin_id) {
                lhs <- h5read(file,  paste0(group_id,"/", i))
                table <- data.frame(data = t(lhs), cell_type = rhs)
                distance_matrix <- vegdist(t(lhs), method = j)
                # distance_matrix <- vegdist(t(lhs), method = "manhattan")
                

                class_cell=data.frame(class=factor(rhs))
                rownames(class_cell)<-rownames(distance_matrix)
                bk<-c(seq(-1,1,by=0.1))
                res_path=paste0('/home/python/higashi/notebook/heatmap/',j,'_',i,'.png')
                # pdf(pdf_path)
                pheatmap_obj<-pheatmap(distance_matrix,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        display_numbers = FALSE,
                        annotation_row = class_cell,
                        annotation_col=class_cell,
                        border_color = 'black',
                        color = colorRampPalette(c('blue', 'white', 'red'))(length(bk)),
                        silent = TRUE  # 防止在图形设备不存在时产生警告
                )
                # dev.off() 
                # 将绘图结果保存为图片
                ggsave(res_path, plot = pheatmap_obj$gtable, device = "png")
                print(res_path)
        }
}

