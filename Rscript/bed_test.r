library(data.table)
lines <- c(
  "# locus chr19:49302001-49304701",
  "# refGene encodeRegions",
  "# zero-based, half-open coords",
  'track name="-log10(shap_Padjust)" description="BedGraph format" visibility=full color=0,153,0 priority=20 plotType="points"'
)
lines1 <- c(
  "# locus chr19:49302001-49304701",
  "# refGene encodeRegions",
  "# zero-based, half-open coords",
  'track name="shap_dist" description="BedGraph format" visibility=full color=0,153,0 priority=20 plotType="points"'
)
bed_path<-"/home/python/higashi/compare_dchic/bed_file/padjust_prebulk_cortex_250k.bedGraph"
dist_path<-"/home/python/higashi/compare_dchic/bed_file/dist_prebulk_cortex_250k.bedGraph"
# 将特定文本写入文件
writeLines(lines, bed_path)
writeLines(lines1, dist_path)
# data<-fread("/home/python/higashi/notebook/cortex250k_compartment_shap_chord_dmanova_w2_default_b150.csv")
# data<-fread("/home/python/higashi/notebook/cortex250k_compartment_shap_chord_dmanova_w2_default_b150_remove_chrX_ncb.csv")
# data<-fread('cortex250k_compartment_shap_chord_dmanova_w2_default_b150_remove_chrX_ncb.csv')


# data<-fread('/home/compartment_raw_sim_shap_manhattan.csv')



# data<-fread('shap_hic.csv')
# data<-fread('/home/compartment_zscore_shap_hic.csv')
data<-fread('/home/python/higashi/notebook/Rscript/compartment_raw_ori_pc1_shap_manhattan.csv')

label<-fread('/home/python/higashi/notebook/Rscript/cortex250k/compartment_raw_bed.csv')


str(label)
filtered_data <- subset(data, effect_size>0.06)
# filtered_data <- subset(data, p_value<10^(-6))

data1=data.frame(matrix(nrow = length(filtered_data$effect_size)))
data1 <- data.frame(p_value = -log10(filtered_data$p_adjust+10^(-200)), index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)
# data1 <- data.frame(p_value = -log10(filtered_data$p_value), index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)

str(data1)

merged_data <- merge(label, data1, by.x = "index", by.y = "index")
data2 <- data.frame(
  chrom = merged_data$chrom,
  start = merged_data$start,
  end = merged_data$end,
  p_value = merged_data$p_value
)

str(data2)
data3=data.frame(matrix(nrow = length(filtered_data$effect_size)))
data3 <- data.frame(dist = filtered_data$F, index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)

merged_data2<-merge(label, data3, by.x = "index", by.y = "index")
data4<- data.frame(
  chrom = merged_data2$chrom,
  start = merged_data2$start,
  end = merged_data2$end,
  dist = merged_data2$dist
)
str(data4)
write.table(data2, bed_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE,append = TRUE)

write.table(data4, dist_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE,append = TRUE)





