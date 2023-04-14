library(metaX)
para <- new("metaXpara")
setwd("/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats")
pfile <- "/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats/cross_species_all_QC_pqn_normalized_QCRSC_pclean_cvfilter_intensity.txt"
sfile <- "/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats/full_sample_list_formetax.txt"
rawPeaks(para) <- read.delim(pfile,check.names = FALSE, sep = '\t')
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
para <- missingValueImpute(para, valueID = "value", method = "knn", cpu = 8)
tab <- getPeaksTable(para, valueID = 'value')

# reformat the table
# reformat normalized data
samples <- as.vector(tab[,1])
name <- colnames(tab)[-1:-4]
tab1 <- tab[,-1:-4]
tab2 <- t(tab1)
tab2 <- as.data.frame(tab2)
tab2 <- cbind(name, tab2)
colnames(tab2) <- c('name', samples)
rownames(tab2) <- 1:dim(tab2)[1]
write.table(tab2, file = '/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats/cross_species_all_QC_pqn_normalized_QCRSC_pclean_cvfilter_KNN_intensity.txt', quote = FALSE, sep = '\t', row.names = FALSE)

