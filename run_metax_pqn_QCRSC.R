library(metaX)
library(stringr)
# run pqn + QC-RSC and write normalized data into a table
para <- new("metaXpara")
setwd("/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats")
pfile <- "/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats/cross_species_filtered_intensity.csv"
sfile <- "/lustre/fengsiyuan/data/metabolome/test_adipose/analysis/personalized_analysis/cross_species_new_data/intensity5000/pos/stats/cross_species_filtered_list.txt"
rawPeaks(para) <- read.delim(pfile,check.names = FALSE, sep = ',')
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
para <- metaX::normalize(para,method="pqn",valueID="value", cpu = 8)
res <- doQCRLSC(para,cpu=8, impute = FALSE)
tab <- getPeaksTable(res$metaXpara, valueID = 'valueNorm')
write.table(tab, file = 'cross_species_pqn_normalized_QCRSC_intensity.txt', quote = FALSE, sep = '\t')

# reformat normalized data
samples <- as.vector(tab[,1])
name <- colnames(tab)[-1:-4]
tab1 <- tab[,-1:-4]
tab2 <- t(tab1)
tab2 <- as.data.frame(tab2)
tab2 <- cbind(name, tab2)
colnames(tab2) <- c('name', samples)
rownames(tab2) <- 1:dim(tab2)[1]

# run pclean and write cleaned data into a table
para <- new("metaXpara")
rawPeaks(para) <- as.data.frame(tab2)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
#dim(para@peaksData)
# para <- dataClean(para, cpu = 8)
tab3 <- getPeaksTable(para)

# feature filter according to the CV (30%) of features in QC sample
samples <- as.vector(tab3[,1])
name <- colnames(tab3)[-1:-4]
tab4 <- tab3[,-1:-4]
tab4 <- t(tab4)
tab4 <- as.data.frame(tab4)
tab4 <- cbind(name, tab4)
colnames(tab4) <- c('name', samples)
rownames(tab4) <- 1:dim(tab4)[1]
tab4[is.na(tab4)] <- 0

qc_index <- which(samples == samples[str_detect(samples,'QC')])
tab5 <- tab4[which(apply(tab4, 1, function(x)all(sd(as.numeric(x[-1][qc_index]))/mean(as.numeric(x[-1][qc_index])) <= 0.3))),]
write.table(tab5, file = 'cross_species_pqn_normalized_QCRSC_pclean_cvfilter_intensity.txt', quote = FALSE, sep = '\t', row.names = FALSE)

