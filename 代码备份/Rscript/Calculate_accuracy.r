library(data.table)
A<-fread('/home/sim_data/meanA-1.csv')
B<-fread('/home/sim_data/meanB-1.csv')
# A<-fread('/home/sim_data/var0.3.csv')
# B<-fread('/home/sim_data/var1.5.csv')
# data<-fread('/home/compartment_raw_sim_var_shap_manhattan.csv')
# data<-fread('compartment_raw_sim_SAWB_shap_manhattan.csv')
data<-fread('/home/compartment_raw_sim_shap_manhattan.csv')
# C<-fread('/home/sim_data/var1.5.csv')
# A$datalabel
# search_term <- "wd"
search_term <- "wd|vd"
resultA <- grep(search_term, A$datalabel)
# resultA 
resultB <- grep(search_term, B$datalabel)

# search_nd <- "nd"
# nA <- grep(search_nd, A$datalabel)
# # resultA 
# nB <- grep(search_nd, B$datalabel)

# combined_nd<-union(nA, nB) #实际未加差异

# resultB
# resultC <- grep(search_term, C$datalabel)
# resultC
# combined_result <- union(union(resultA, resultB),resultC)

combined_result <- union(resultA, resultB) #实际加了差异


N<-length(data$effect_size)
print(paste("N",N))
# data<-fread('/home/compartment_raw_sim_var_shap_manhattan.csv')
thresholds<-0.16
print(paste("thresholds",thresholds))
filtered_data <- subset(data, effect_size>thresholds)
#提取bin编号
filtered_data$bin_number <- as.numeric(sub("bin_", "", filtered_data$bin))#预测为有显著差异
data$bin_number <- as.numeric(sub("bin_", "", data$bin))
combined_nd<-setdiff(data$bin_number,combined_result) 
N_number<-setdiff(data$bin_number,filtered_data$bin_number) #预测为无显著差异
length(N_number)
length(filtered_data$bin_number)
length(combined_nd)
length(combined_result)
# 合并resultA和resultB并保留唯一值
# combined_result <- union(resultA, resultB)
TP=length(intersect(filtered_data$bin_number, combined_result))
TN=length(intersect(combined_nd, N_number))
FN=length(intersect(combined_result, N_number))
FP=length(intersect(filtered_data$bin_number,combined_nd))

# accuracy <- length(intersect(filtered_data$bin_number, combined_result)) / length(combined_result)

print(paste("TP:", TP))
print(paste("TN:", TN))
print(paste("FP:",FP))
print(paste("FN:", FN))
precision=TP/(TP+FP)
recall=TP/(TP+FN)
acc=(TP+TN)/N
f1=(2*precision*recall)/(precision+recall)
FPR=FP/(TN+FP)
print(paste0("precision:",format(precision*100,digits = 3),"%"))
# print(paste("precision:", precision))
print(paste0("recall:", format(recall*100,digits = 3),"%"))
# print(paste("recall:", recall))
# print(paste("FPR:", FPR))
print(paste0("FPR:", format(FPR*100,digits = 3),"%"))
# print(paste("acc:",acc))
print(paste0("acc:",format(acc*100,digits = 3),"%"))
print(paste("f1:", f1))




# library(pROC)

# library(data.table)
# A<-fread('/home/sim_data/meanA-1.csv')
# B<-fread('/home/sim_data/meanB-1.csv')
# # A<-fread('/home/sim_data/var0.3.csv')
# # B<-fread('/home/sim_data/var1.5.csv')
# # data<-fread('/home/compartment_raw_sim_var_shap_manhattan.csv')
# # data<-fread('compartment_raw_sim_SAWB_shap_manhattan.csv')
# data<-fread('/home/compartment_raw_sim_shap_manhattan.csv')
# data$bin_number <- as.numeric(sub("bin_", "", data$bin))
# # C<-fread('/home/sim_data/var1.5.csv')
# # A$datalabel
# # search_term <- "wd"
# search_term <- "wd|vd"
# resultA <- grep(search_term, A$datalabel)
# # resultA 
# resultB <- grep(search_term, B$datalabel)

# # search_nd <- "nd"
# # nA <- grep(search_nd, A$datalabel)
# # # resultA 
# # nB <- grep(search_nd, B$datalabel)

# # combined_nd<-union(nA, nB) #实际未加差异

# # resultB
# # resultC <- grep(search_term, C$datalabel)
# # resultC
# # combined_result <- union(union(resultA, resultB),resultC)

# combined_result <- union(resultA, resultB) #实际加了差异
# combined_nd<-setdiff(data$bin_number,combined_result)


# # 示例数据（替换为你自己的数据）
# predicted_scores <- data$effect_size
# # true_labels <- ifelse(combined_nd %in% data$bin_number , 1, 0)
# true_labels <- as.integer(data$bin_number %in% combined_nd)

# max_td<-max(data$effect_size)

# # 设置不同的阈值
# # thresholds <- seq(0.16, 0.06, by = -0.01)
# thresholds <- seq(0.16, max_td, by = 0.01)
# # 初始化存储真正例率和假正例率的向量
# tpr <- numeric(length = length(thresholds))
# fpr <- numeric(length = length(thresholds))

# # 计算不同阈值下的真正例率和假正例率
# for (i in seq_along(thresholds)) {

#     threshold <- thresholds[i]

#     # 预测标签
#     predicted_labels <- ifelse(predicted_scores > threshold, 1, 0)

#     #   # 计算混淆矩阵
#     #   confusion_matrix <- table(predicted_labels, true_labels)

#     # 计算真正例率和假正例率
#     # 计算同时为1的数量
#     both_ones <- sum(predicted_labels == 1 & true_labels == 1)

#     # 计算预测为1但真实为0的数量
#     predicted_one_actual_zero <- sum(predicted_labels == 1 & true_labels == 0)
#     tpr[i] <- both_ones/sum(true_labels == 1)
#     fpr[i] <- predicted_one_actual_zero/sum(true_labels == 0)
#     # roc_curve <- roc(true_labels, predicted_labels)


# }


# pdf('AUC.pdf')
# # 绘制 ROC 曲线
# plot(fpr,tpr, main = "ROC Curve", col = "blue", lwd = 2)

# # # 添加 AUC 值的图例
# # legend("bottomright", legend = paste("AUC =", round(roc_auc, 2)), col = "blue", lwd = 2)
# dev.off()