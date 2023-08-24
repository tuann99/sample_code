# AS.410.671.81 – Final Project
# Tuan Nguyen

# *** Primary Reference *** : https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3227

# Title:
#	- Social isolation stress-induced aggression model
# Summary:
#	- Analysis of heads of Canton-S males subjected to social isolation (single-housed).
# 	  Social experience with conspecifics is an environmental influence on aggressiveness common to many species.
#     Results provide insight into the molecular basis of the effect of social experience on aggressiveness.
# Organism:
#	- Drosophila melanogaster
# Platform:
#	- 	GPL1322: [Drosophila_2] Affymetrix Drosophila Genome 2.0 Array
#	- 	Annotation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1322
# Citation:
#	- Wang et al. "A common genetic target for environmental and heritable influences on aggressiveness in Drosophila."
#	- Proc Natl Acad Sci U S A 2008 Apr 15;105(15):5657-63.
#	- PMID: 18408154
#	- Reference Series: GSE6994

# set working directory
wd <- "C:\\Users\\tuann\\OneDrive - Johns Hopkins\\School\\AS.410.671.82 (Gene Expression Data Analysis and Visualization)\\final_project"
setwd(wd)

# installations if needed. Some packages are old and require manual installation of package and dependencies.
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GEOquery")
# install.packages("stats")
# install.packages("outliers")
# untar("EMA_1.4.7.tar.gz")
# untar("heatmap.plus_1.3.tar.gz")
# BiocManager::install("siggenes", force = TRUE)
# BiocManager::install("biomaRt")
# install.packages("GSA")
# install.packages("FactoMineR")
# install.packages("heatmap.plus", repos = NULL, type = "source")
# install.packages("EMA", repos = NULL, type = "source")
# install.packages("e1071")
# install.packages("scatterplot3d")

# load libraries
library(GEOquery)
library(stats)
library(outliers)
library(EMA)
library(e1071)
library(MASS)
library(multtest)
library(kernlab)
library(scatterplot3d)
library(R.utils)

# Must load SVM-RFE.r from source to perform SVM feature ranking.
# Source can be downloaded to local machine via:
# https://github.com/johncolby/SVM-RFE/blob/master/msvmRFE.R
source(file = "msvmRFE.R")

# Reading in GEO Data
gds <- getGEO("GDS3227")
exp_data <- Table(gds)

# Reading in the corresponding GEO annotation file
# ann_file <- getGEO("GPL1322", AnnotGPL = TRUE) # didn't work so I downloaded the file manually
# gunzip("GPL1322.annot.gz", remove=FALSE)
ann_file <- "GPL1322.annot"
ann_data <- read.delim(ann_file, header = TRUE, row.names = 1, skip = 27, sep = "\t")
ann_data <- ann_data[1:nrow(exp_data), ]
gene_info <- data.frame(Description = ann_data$Gene.title, Symbol = ann_data$Gene.symbol)
rownames(gene_info) <- rownames(ann_data)

# Appending Group (Grp) and Single Housing (Sng) to corresponding samples
for (i in 3:5) {
	names(exp_data)[i] <- paste("Sng", names(exp_data[i]), sep = "_")
}
for (j in 6:8) {
	names(exp_data)[j] <- paste("Grp", names(exp_data[j]), sep = "_")
}
# colnames(exp_data) # check names
samp_matrix <- data.matrix(exp_data[, (3:ncol(exp_data))])
rownames(samp_matrix) <- rownames(ann_data)
samp_count <- ncol(samp_matrix)
profile_count <- nrow(samp_matrix)

# Descriptive Satistics
data_sd <- apply(samp_matrix, 1, sd, na.rm = TRUE)
data_rowMeans <- rowMeans(samp_matrix, na.rm = TRUE)

# Histograms illustrating initial spread
png("Histogram_data_rowMeans.png")
hist(
	data_rowMeans, 
	col  = "Red",
	xlab = "Mean expression value for Sng/Grp samples",
	ylab = "Frequency",
	main = paste("Histogram of mean expression values for", profile_count, "profiles")
	)
dev.off()

png("Histogram_data_stdev.png")
hist(
	data_sd,
	col  = "Blue",
	xlab = "Standard Deviation expression value for Sng/Grp samples",
	ylab = "Frequency",
	main = paste("Histogram of standard deviation expression values for", profile_count, "profiles")
	)
dev.off()

########################
#   Outlier Analysis   #
########################

# Correlation Matrix
cor_matrix  <- cor(samp_matrix, method = "pearson", use = "pairwise.complete.obs")

color <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000",
		   "#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")

png("Heatmap_cor_matrix.png", width = 15, height = 15, units = "in", res = 300)
par(mar = c(10, 3, 3, 10))
heatmap(
	cor_matrix,
	col = color,
	scale = "column",
	xlab = "Sample Data",
	ylab="Sample Data",
	main = "Heatmap illustrating expression values of Sng/Grp samples"
	)
dev.off()

# Column Means v. Column Variances
col_mean <- apply(log2(samp_matrix), 2, mean)
col_var <- apply(log2(samp_matrix), 2, var)
cv <- col_var / col_mean
png("Scatterplot_ColMeansCV.png")
plot(
	col_mean,
	cv,
	xlab = "log2(ColMean)",
	ylab = "log2(CV)",
	main = "Plot of Column Mean v. Column Variance for Sng/Grp samples",
	col  = c(rep("Green", samp_count/2), rep("Blue", samp_count/2)),
	pch  = c(rep(17, samp_count/2), rep(19, samp_count/2))
	)
legend("topright", c("Single Housed", "Group Housed"), pch = c(17, 19), col = c("Green", "Blue"))
text(col_mean, cv, labels = names(col_mean), cex = 0.5, offset = 10)
dev.off()

# Row Means v. Column Variances
row_mean <- apply(log2(samp_matrix), 1, mean)
row_var  <- apply(log2(samp_matrix), 1, var)
r_cv     <- row_var / row_mean
# head(r_cv)
png("Scatterplot_RowMeansCV.png")
plot(
	row_mean,
	r_cv,
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Plot of Row Mean v. Row Variance for Sng/Grp samples",
	col  = c(rep("Green", samp_count/2), rep("Blue", samp_count/2)),
	pch  = c(rep(17, samp_count/2), rep(19, samp_count/2))
	)
legend("topright", c("Single Housed", "Group Housed"), pch = c(17, 19), col = c("Green", "Blue"))
dev.off()

# Correlation plot of Row Averages
cor_means <- apply(cor_matrix, 1, mean)

png("Scatterplot_CorMeans.png")
plot(
	c(1,length(cor_means)),
	range(cor_means),
	type = "n",
	xlab = "",
	ylab = "Average correlation",
	main = "Avg correlation for Sng/Grp samples",
	axes = FALSE
	)
points(
	cor_means,
	col = c(rep("Green", samp_count/2), rep("Blue", samp_count/2)),
	pch = c(rep(17, samp_count/2), rep(19, samp_count/2))
	)
axis(1, at=c(1:length(cor_means)), labels = colnames(samp_matrix), las = 2, cex.lab = 0.4, cex.axis = 0.6)
axis(2)
grid(nx = 16, col = "grey")
legend(
	"topright", 
	c("Single Housed", "Group Housed"), 
	pch = c(17, 19), col = c("Green", "Blue"), bg = "white"
	)
dev.off()

# Identifying Outlier(s) via outlier()
o <- cor_means <= outlier(cor_means)
outlier <- cor_means[o]
cat(sprintf("%s Outlier(s) identified\n", length(outlier)))
# outlier
# 1 Outlier(s) identified
# > outlier
# Grp_GSM161262
#     0.9844887

# Remove Outlier(s)
data_no_outliers <- rm.outlier(cor_means, fill = FALSE, median = FALSE, opposite = FALSE)

########################
#     Filter Genes     #
########################
quantile(log2(rowMeans(samp_matrix)))
# Output:
#        0%       25%       50%       75%      100%
# -1.454131  3.088798  4.825013  7.425855 13.303082

# Eliminating probes with rowMeans less than 0 on a log2 scale
samp_matrix_filtered <- subset(samp_matrix, log2(rowMeans(samp_matrix)) > 0)
removed <- nrow(samp_matrix) - nrow(samp_matrix_filtered)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed)) # 55

# Use expFilter() to filter genes with low expression values
# A gene is kept if at least 0.01*ncol(samp_matrix) of its values is higher than threshold (3).
dat_fil <- expFilter(log2(samp_matrix_filtered), threshold = 3, p = 0.01, graph = TRUE)
dat_fil <- subset(dat_fil, rowMeans(dat_fil) > 0)
num_lowexp <- nrow(samp_matrix_filtered) - nrow(dat_fil)
cat(sprintf("%s gene(s) identified and removed for low expression\n", num_lowexp)) # 1308

# Row Means v. Column Variances on filtered data
fil_mean <- apply(dat_fil, 1, mean)
fil_var <- apply(dat_fil, 1, var)
f_cv <- fil_var / fil_mean

# Plotting filtered genes
png("Scatterplot_RowMeansCV_filt_Stage1.png")
plot(
	fil_mean,
	f_cv,
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Plot of Row Mean v. Row Variance for Filtered Sng/Grp samples",
	col  = c(rep("Green", samp_count/2), rep("Blue", samp_count/2)),
	pch  = c(rep(17, samp_count/2), rep(19, samp_count/2))
	)
legend("topright", c("Single Housed", "Group Housed"), pch = c(17, 19), col = c("Green", "Blue"))
abline(v = 3, col = 2, lwd = 2) # Threshold determined for stage 2/2 of the filtering process
dev.off()

# Eliminating probes with rowMeans less than 3 on a log2 scale
dat_filtered <- subset(dat_fil, rowMeans(dat_fil) > 3)
removed_2 <- nrow(dat_fil) - nrow(dat_filtered)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed_2)) # 4017
# nrow(dat_filtered) # 13572

fil_mean_2 <- apply(dat_filtered, 1, mean)
fil_var_2 <- apply(dat_filtered, 1, var)
fil_cv_2 <- fil_var_2 / fil_mean_2

png("Scatterplot_RowMeansCV_filt_Stage2.png")
plot(
	fil_mean_2,
	fil_cv_2,
	xlim = c(0, 12), ylim = c(0, 12), # Consistant window Size w/ previous plot
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Stage 2/2 Filtered Plot of Row Mean v. Row Variance for Sng/Grp samples",
	col  = c(rep("Green", samp_count/2), rep("Blue", samp_count/2)),
	pch  = c(rep(17, samp_count/2), rep(19, samp_count/2))
	)
legend("topright", c("Single Housed", "Group Housed"), pch = c(17, 19), col = c("Green", "Blue"))
abline(v = 3, col = 2, lwd = 2)
dev.off()

# Update gene_info dataframe
gene_info <- subset(gene_info, rownames(gene_info) %in% rownames(dat_filtered))

########################
#  Feature Selection   #
########################   

type <- lapply(colnames(dat_filtered),
	function(x) {
		if(regexpr("Sng", x) < 1) {"Grp"} else {"Sng"}
	})

# Determine individual P-values in the original expression scale (not log2)
pval <- c()
ptm <- proc.time()
for (i in seq(nrow(dat_filtered))){
	t <- t.test(2^dat_filtered[i, type == "Sng"], 2^dat_filtered[i, type == "Grp"], alternative = "two.sided")
	pval <- c(pval, t$p.value)
}
proc.time() - ptm
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(pval < 0.05))) # 581
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(pval < 0.01))) # 80

Bonferroni <- 0.05/length(pval)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni)) # 0
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(pval < Bonferroni))) # 0

# Adjust P-values via Holm's method
p_holm <- p.adjust(pval, method = "holm")
p_raw  <- sort(pval)
p_holm <- sort(p_holm)
p1 <- as.matrix(p_raw)
p2 <- as.matrix(p_holm)
all_p <- cbind(p1, p2)
colnames(all_p) <- c("Raw P-value", "Adjusted P-Value")

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(p_holm < 0.05))) # 0
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(p_holm < 0.01))) # 0

Bonferroni_2 <- 0.05/length(p_holm)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni_2))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(pval < Bonferroni_2))) # 0

# Plot the sorted raw P-values
png("Raw_PvaluePlot.png")
plot(
	p_raw, type = "b", pch = 1, col = "Pink",
	xlab = "Genes",
	ylab = "P-values",
	main = "Raw P-values for All Genes\n"
	)
dev.off()
	
# Plot Adjusted vs. Raw P-values
png("PvaluePlot.png")
matplot(
	all_p, type = "b", pch = 1, col = 1:2,
	xlab = "Genes",
	ylab = "P-values",
	main = "Adjusted Vs. Raw P-values for All Genes\n"
	)
legend("bottomright", legend = colnames(all_p), pch = 1, col = 1:2)
dev.off()

gene_info$pvalue <- pval

# Hypothesis Pass/Fail Dataframe
thresh <- 0.05
rnames <- rownames(dat_filtered)
p_test <- lapply(as.list(pval), function(f) { 
	if (f < thresh) TRUE else FALSE })
pval_df <- as.data.frame(do.call(rbind, p_test), rname = rnames)
names(pval_df) <- c("Pval_Test")
rownames(pval_df) <- rownames(dat_filtered)
gene_info$Pval_Test <- pval_df$Pval_Test

cat(sprintf("Threshold: %s\n", thresh))
table(pval_df$Pval_Test)
#   Number of True        : 581
#   Number of False       : 12991
#   Total (nrow(pval_df)) : 13572

# Histogram illustrating the distribution of P-values. Q4
png("Histogram_pval.png")
hist(
	pval,
	col  = "Pink",
	xlab = "P-Value",
	ylab = "Frequency",
	main = "Histogram of T-Test P-Values for GDS3227 profiles"
	)
dev.off()

# Determining fold change between treatments in log2
Sng_matrix <- dat_filtered[, type == "Sng"]
Sng_mean <- apply(Sng_matrix, 1, mean, na.rm = TRUE)
Grp_matrix <- dat_filtered[, type == "Grp"]
Grp_mean <- apply(Grp_matrix, 1, mean, na.rm = TRUE)

fold <- Sng_mean - Grp_mean
# 2^max(fold) # 97.54693
# 2^min(fold) # 0.009346909

# Update gene_info data_frame with fold information
gene_info$fold <- fold

# probesets that demonstrate a 2x fold change
fold_test <- lapply(as.list(fold),
	function(x) {
		if (abs(x) > log2(2)) TRUE else FALSE })
fold_df <- as.data.frame(do.call(rbind, fold_test))
names(fold_df) <- c("Fold_Test")

# Update gene_info data_frame with fold information
gene_info$Fold_Test <- fold_df$Fold_Test

table(fold_df$Fold_Test)
# Number of True        : 481
# Number of False       : 13091
# Total (nrow(fold_df)) : 13572

# Histogram illustrating the distribution of log2(fold change)
png("Histogram_fold.png")
hist(
	fold,
	col  = "Lightblue",
	xlab = "Log2 Fold Change",
	ylab = "Frequency",
	main = paste("Histogram of Fold Change values for GDS3227 profiles")
	)
abline(v = log2(2), col = 2, lwd = 2)
abline(v = -log2(2), col = 2, lwd = 2)
dev.off()

# overall Volcano plot showing cutoffs and differences within the dataset
p_transformed <- (-1 * log10(pval))
png("VolcanoPlot.png")
plot(
	range(p_transformed),
	range(fold),
	type = "n", xlab = "-1 * log10(P-Value)", ylab = "Fold Change",
	main = "Volcano Plot From Single Housed and Group Housed Drosophila"
	)
points(
	p_transformed, fold,
	col = 1, bg = 1, pch = 21
	)
points(
	p_transformed[(p_transformed > -log10(0.05) & fold > log2(2))],
	fold[(p_transformed > -log10(0.05) & fold > log2(2))],
	col = 1, bg = 2, pch = 21
	)
points(
	p_transformed[(p_transformed > -log10(0.05) & fold < -log2(2))],
	fold[(p_transformed > -log10(0.05) & fold < -log2(2))],
	col = 1, bg = 3, pch = 21
	)
abline(v = -log10(0.05))
abline(h = -log2(2))
abline(h = log2(2))
dev.off()

passed_genes <- subset(gene_info, (Fold_Test & Pval_Test) == TRUE)
cat(sprintf("Total number of genes that pass both (Pval and Fold) tests: %s\n", nrow(passed_genes))) # 76

# write genes with their corresponding values to a file
# write.table(passed_genes, file = "passedgenes.csv", sep = ",", col.names = NA, qmethod = "double")

# Ordering the highest genes (by P-value)
# Note: dat_filtered is still in log2 scale
best_genes <- order(pval)[1:length(pval)]
best_genes_df <- data.frame(
	index = best_genes, exp = 2^dat_filtered[best_genes, ],
	pval = pval[best_genes])

# Expression matrix with the 'best genes' in the original scale (based on P-value ranking)
top_genes_matrix <- 2^dat_filtered[best_genes, ]

# Feature Selection via svmRFE which utilizes the library e1071
t_dat <- t(top_genes_matrix) # transpose to get the correct format
label <- as.vector(unlist(type))
svm_df <- data.frame(label, t_dat)
svm_df$label <- as.factor(ifelse(svm_df$label == "Sng", 0, 1)) # 0 = Sng, 1 = Grp
ranked_list <- svmRFE(svm_df, k = 3, halve.above = 1000)
top_ranked_genes <- top_genes_matrix[ranked_list, ]
rownames(top_ranked_genes) <- rownames(top_genes_matrix[ranked_list, ])

# Create a new gene info dataframe for the ranked genes
top_genes_info <- gene_info[rownames(top_ranked_genes), ]
tg <- top_genes_info$pvalue[top_genes_info$pvalue < thresh]

png("Histogram_ranked_pval.png")
hist(
	tg,
	col  = "Pink",
	xlab = "P-Value",
	ylab = "Frequency",
	main = "Histogram of Ranked T-Test P-Values"
	)
abline(v = thresh, col = 2, lwd = 2)
dev.off()

top5 <- head(top_genes_info, n = 5L)
bottom5 <- tail(top_genes_info, n = 5L)
# top5
# bottom5
# > top5
#                                                Description Symbol       pvalue
# 1637077_s_at CG3777 gene product from transcript CG3777-RB CG3777 0.0008974962
# 1631441_s_at                                      slowdown   slow 0.0005429153
# 1637532_s_at                 Dipeptidyl aminopeptidase III DppIII 0.0035840084
# 1640155_at                     Inositol-requiring enzyme-1   Ire1 0.0009024727
# 1624121_at                                       knickkopf    knk 0.0004706164

# > bottom5
#                                                  Description  Symbol    pvalue
# 1628486_a_at   CG2682 gene product from transcript CG2682-RB      d4 0.8307612
# 1629081_at                             orientation disruptor     ord 0.5232780
# 1631323_a_at   CG2663 gene product from transcript CG2663-RA  CG2663 0.6576469
# 1639794_at   CG31812 gene product from transcript CG31812-RB CG31812 0.7191000
# 1636451_a_at               Activating transcription factor-2   Atf-2 0.8447166

###########################
#  Sample Classification  #
###########################  

# PCA analysis on the first 3 component vectors
pca <- prcomp(top_ranked_genes)
pca_loads <- pca$x[, 1:3]
label <- as.factor(ifelse(label == "Sng", 0, 1)) # 0 = Sng, 1 = Grp
svp <- ksvm(pca_loads, label, type = "C-svc")
# svp
# Support Vector Machine object of class "ksvm"

# SV type: C-svc  (classification)
#  parameter : cost C = 1

# Gaussian Radial Basis kernel function.
#  Hyperparameter : sigma =  55.8161537365673

# Number of Support Vectors : 13136

# Objective Function Value : -11852.46
# Training error : 0.382847

# Get fitted values
fit <- fitted(svp)

# error rates (incorrect classifications)
er1 <- sum(fit[label == "Sng"] == "Grp") # 0
er2 <- sum(fit[label == "Grp"] == "Sng") # 0

# plot the data
png("Scatterplot_PCA_1v2.png")
plot(
	range(pca_loads[, 1]),
	range(pca_loads[, 2]),
	type = "n",
	xlab = "Principal Component 1",
	ylab = "Principal Component 2",
	main = "PCA Plot for GDS3227 data\n PC1 vs. PC2"
	)
points(
	pca_loads[, 1][as.numeric(type == "Sng") == 1],
	pca_loads[, 2][as.numeric(type == "Sng") == 1],
	col = "Red", pch = 15
	)
points(
	pca_loads[, 1][as.numeric(type == "Grp") == 1],
	pca_loads[, 2][as.numeric(type == "Grp") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomleft", 
	c("Single Housed", "Group Housed"),
	col = c("Blue", "Red"), pch = c(19, 15)
	)
dev.off()
	
png("Scatterplot_PCA_1v3.png")
plot(
	range(pca_loads[, 1]), 
	range(pca_loads[, 3]), 
	type = "n",
	xlab = "Principal Component 1",
	ylab = "Principal Component 3",
	main = "PCA Plot for GDS3227 data\n PC1 vs. PC3"
	)
points(
	pca_loads[, 1][as.numeric(type == "Sng") == 1],
	pca_loads[, 3][as.numeric(type == "Sng") == 1],
	col = "Red", pch = 15
	)
points(
	pca_loads[, 1][as.numeric(type == "Grp") == 1],
	pca_loads[, 3][as.numeric(type == "Grp") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomright",
	c("Single Housed", "Group Housed"),
	col = c("Blue", "Red"), pch = c(19, 15)
	)
dev.off()

png("Scatterplot_PCA_2v3.png")
plot(
	range(pca_loads[, 2]), 
	range(pca_loads[, 3]), 
	type = "n",
	xlab = "Principal Component 2",
	ylab = "Principal Component 3",
	main = "PCA Plot for GDS3227 data\n PC2 vs. PC3"
	)
points(
	pca_loads[, 2][as.numeric(type == "Sng") == 1], 
	pca_loads[, 3][as.numeric(type == "Sng") == 1],
	col = "Red", pch = 15
	)
points(
	pca_loads[, 2][as.numeric(type == "Grp") == 1], 
	pca_loads[, 3][as.numeric(type == "Grp") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomleft", 
	c("Single Housed", "Group Housed"), 
	col = c("Blue", "Red"), pch = c(19,15)
	)
dev.off()
	
pca_var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
cat(sprintf("Note: Approximately %s variability is explained using the only the first two eigenvalues.\n", 
	sum(pca_var[1:2])))

# Scree Plot adapted from Lab_07 Illistrating Level of Variance
png("Screeplot.png")
plot(
	c(1:length(pca_var)), 
	pca_var, 
	type = "b", 
	xlab = "Components",
	ylab = "Percent Variance", 
	bg = "Blue", pch = 21
	)
title("Scree Plot Illistrating %-Variability Explained By Each Eigenvalue")
dev.off()

# MDS Analysis via Kruskal’s Non-metric Approach
dat_dist <- dist(t(top_ranked_genes))
dat_mds <- isoMDS(dat_dist)

png("MDSplot_Kruskal.png")
plot(dat_mds$points, type = "n")
points(
	dat_mds$points[, 1][as.numeric(type == "Sng") == 1],
	dat_mds$points[, 2][as.numeric(type == "Sng") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	dat_mds$points[, 1][as.numeric(type == "Grp") == 1], 
	dat_mds$points[, 2][as.numeric(type == "Grp") == 1], 
	col = "Blue", pch = 16, cex = 1.5)
title(main = "Kruskal’s Non-metric MDS plot for GDS3227 Dataset\nSingle Housed vs. Group Housed")
legend("bottomleft", c("Single Housed", "Group Housed"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

# MDS Analysis via the Classical Metric Approach
dat_loc <- cmdscale(dat_dist)

png("MDSplot_Classical.png")
plot(dat_loc, type = "n")
points(
	dat_loc[, 1][as.numeric(type == "Sng") == 1],
	dat_loc[, 2][as.numeric(type == "Sng") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	dat_loc[, 1][as.numeric(type == "Grp") == 1], 
	dat_loc[, 2][as.numeric(type == "Grp") == 1], 
	col = "Blue", pch = 16, cex = 1.5)
title(main = "Classical Metric MDS plot for GDS3227 Dataset\nSingle Housed vs. Group Housed")
legend("bottomleft", c("Single Housed", "Group Housed"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

# Weighted Laplacian Graph
# Note: The function GeneratePhi was taken from Lab_07 and Lecture_08
GeneratePhi <- function (X, qnt = NULL) {
	dist2full <- function(dis) {
		n <- attr(dis, "Size")
		full <- matrix(0, n, n)
		full[lower.tri(full)] <- dis
		full + t(full)
	}
	dat.dis <- dist(t(X),"euc")^2
	if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
	if(is.null(qnt))  {eps <- min(dat.dis[dat.dis != 0])}
	kernal <- exp(-1 * dat.dis / (eps))
	K1 <- dist2full(kernal)
	diag(K1) <- 0
	D <- matrix(0, ncol = ncol(K1), nrow = ncol(K1))
	tmpe <- apply(K1, 1, sum)
	tmpe[tmpe > 0] <- 1/sqrt(tmpe[tmpe > 0])
	tmpe[tmpe < 0] <- 0
	diag(D) <- tmpe
	L <- D%*% K1 %*% D
	X <- svd(L)$u
	Y <- X / sqrt(apply(X^2, 1, sum))
}

temp <- t(top_ranked_genes)
temp <- scale(temp, center = TRUE, scale = TRUE)
phi <- GeneratePhi(t(temp), qnt = NULL)

png("LaplacianPlot.png")
plot(
	range(phi[, 1]), range(phi[, 2]),
	xlab = "Phi 1", ylab = "Phi 2",
	main = "Weighted Graph Laplacian Plot for GDS3227 Dataset\nSng vs. Grp"
	)
points(
	phi[, 1][as.numeric(type == "Sng") == 1],
	phi[, 2][as.numeric(type == "Sng") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	phi[, 1][as.numeric(type == "Grp") == 1], 
	phi[, 2][as.numeric(type == "Grp") == 1], 
	col = "Blue", pch = 16, cex = 1.5
	)
legend("top", c("Single Housed", "Group Housed"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

########################
#   Cluster Analysis   #
########################

# Hierarchical Clustering via Manhattan 
top_dist <- dist(t(top_ranked_genes), method = "manhattan")
top_clus <- hclust(top_dist, method = "median")

png("Dendogram_TopGenes.png")
plot(
	top_clus,
	labels = colnames(top_ranked_genes),  
	xlab   = "Clustered Samples",
	ylab   = "Distance",
	main   = "Hierarchical Clustering Dendrogram"
	)
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, the top 50 ranked genes were plotted
png("Heatmap_all_top_genes.png")
heatmap(
	top_ranked_genes[1:nrow(top_ranked_genes), ],
	col  = color,
	xlab = "Samples",
	ylab = "Top Ranked Genes",
	main = "Heatmap for all top ranked genes"
	)
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, 50 randomly ranked genes were plotted 
png("Heatmap_TopRandom50Genes.png")
heatmap(
	top_ranked_genes[sample(top_ranked_genes, 50), ],
	col  = color,
	xlab = "Samples",
	ylab = "Randomly Ranked Genes",
	main = "Heatmap for 50 randomly ranked genes"
	)
dev.off()

# K-means clustering via PCA Analysis 
# Done though the kernlab library
# Extract out the first 5 component vectors and compute K-means with two centers 
dat_kpca <- kpca(t(top_ranked_genes), kernel = "rbfdot", kpar = list(sigma = 0.002), features = 5)
pcv <- pcv(dat_kpca)
rot <- rotated(dat_kpca)
pcv_k <- kmeans(pcv, centers = 2, iter.max = 20)
rot_k <- kmeans(rot, centers = 2, iter.max = 20)

# 2D scatterplot of the first 5 eigenfeatures (from PCA)
png("Scatterplot_PCA.png")
par(mfrow = c(2, 1))
plot(
	pcv, col = pcv_k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC vectors (column-wise)\nk = 2",
	xlab = "P1",
	ylab = "P2"
	)
points(pcv_k$centers, col = 1:2, pch = "*", cex = 2.5)

plot(
	rot, col = rot_k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC Projected vectors\nk = 2",
	xlab = "P1",
	ylab = "P2"
	)
points(rot_k$centers, col = 1:2, pch = "*", cex = 2.5)
dev.off()