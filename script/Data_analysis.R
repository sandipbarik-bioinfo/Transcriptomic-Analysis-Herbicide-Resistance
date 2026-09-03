# Installing BiocManager
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Using BiocManager to install required packages
BiocManager::install("affy",force = TRUE)
BiocManager::install("affycoretools",force = TRUE)
BiocManager::install("arrayQualityMetrics",force = TRUE)
BiocManager::install("limma",force = TRUE)
BiocManager::install("org.Hs.eg.db",force = TRUE)
BiocManager::install("clusterProfiler",force = TRUE)
BiocManager::install("ggplot2",force = TRUE)
BiocManager::install("pheatmap",force = TRUE)
BiocManager::install("hgu133plus2.db",force = TRUE)
BiocManager::install("R2HTML",force = TRUE)
BiocManager::install("annotate",force = TRUE)
BiocManager::install("dplyr",force = TRUE)
BiocManager::install("XML",force = TRUE)
BiocManager::install("ggfortify",force = TRUE)
BiocManager::install("openxlsx",force = TRUE)
BiocManager::install("splitstackshape",force = TRUE)
BiocManager::install("GEOquery",force = TRUE)
BiocManager::install("DEGseq",force = TRUE)
BiocManager::install("qvalue",force = TRUE)
BiocManager::install("Biobase",force = TRUE)
BiocManager::install("DESeq2",force = TRUE)
BiocManager::install("tidyverse",force = TRUE)
BiocManager::install("airway",force = TRUE)
BiocManager::install("RColorBrewer",force = TRUE)
BiocManager::install("PoiClaClu",force = TRUE)
BiocManager::install("metaMA",force = TRUE)
BiocManager::install("metaRNASeq",force = TRUE)


# Defining our working directory

setwd("D:\\Bioinformatics\\GEO_Datasets\\sample_data2")

#Reading the text file

raw_data <- read.table("GSE204857_raw_counts.txt")
print(raw_data)

# Convert to a data frame
df_raw <- data.frame(raw_data)

# Download the dataset
gse <- getGEO("GSE204857")

# Extract the phenotype data associated with the GEO dataset.
metadata <- pData(phenoData(gse[[1]]))
metadata

#selected specific columns
selected_columns <- metadata [, c( "title","genotype:ch1", "treatment:ch1")]

#Creating new data frame
new_metadata <- data.frame(selected_columns)

# Write metadata table to a .txt file
ab <- write.table(new_metadata, file = "metadata_main.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Assuming your dataframes are called raw_data and df_raw

row_names <- df_raw[, 1]  # Extract values from the first column of df_raw
row.names(raw_data) <- row_names  # Assign these values as row names to raw_data
# Remove the first column from raw_data
raw_data <- raw_data[, -1]
# Print the structure of the dataframe
str(raw_data)
# Modify the first row name to an empty string
rownames(raw_data)[1]<-""
colnames(raw_data)
# Assign the values of the first row to column names
colnames(raw_data) <- as.character(raw_data[1, ])


# Remove the first row (since it's used for column names now)
raw_data <- raw_data[-1, ]

# Print the structure of the dataframe
str(raw_data)
# Get the sorted indices of column names
sorted_indices <- order(colnames(raw_data))

# Rearrange the columns according to the sorted indices
raw_data <- raw_data[, sorted_indices]

# Print the structure of the dataframe
str(raw_data)


# Get the sorted indices of the first column
sorted_indices <- order(new_metadata[, 1])

# Rearrange the rows according to the sorted indices
new_metadata <- new_metadata[sorted_indices, ]

# Assign the values of the first column to column names
rownames(new_metadata) <- new_metadata[, 1]
rownames(new_metadata)
# Remove the first column (since it's now used for column names)
new_metadata <- new_metadata[, -1]


newmetatxt <- write.table(new_metadata, file = "metadata_final.txt", sep = "\t", quote = FALSE, row.names = TRUE)
rawsort <- write.table(raw_data, file = "raw_data_sorted.txt", sep = "\t", quote = FALSE, row.names = TRUE)

dim(raw_data)
dim(new_metadata)

# making sure the row names in colData matches to column names in counts_data

all(rownames(new_metadata) %in% colnames(raw_data))

# are they in the same order?

all(rownames(new_metadata) == colnames(raw_data))
class(raw_data)
typeof(raw_data)

# Convert data frame to matrix
raw_data_matrix <- data.matrix(raw_data)
class(raw_data_matrix)


# Subset raw_data to include only columns corresponding to control samples

control_raw_data <- raw_data[, grepl("^c", colnames(raw_data))]

# Subset raw_data to include only columns corresponding to treated samples

treated_raw_data <- raw_data[, grepl("^t", colnames(raw_data))]


# Subset raw_data for susceptible samples (both control and treated)
susceptible_cols <- grep("S", colnames(raw_data), value = TRUE)
susceptible_raw_data <- raw_data[, susceptible_cols]

# Subset raw_data for resistant samples (both control and treated)
resistant_cols <- grep("R", colnames(raw_data), value = TRUE)
resistant_raw_data <- raw_data[, resistant_cols]

# Subset new_metadata to include only samples with treatment condition "control"

new_metadata_control <- subset(new_metadata, treatment.ch1 == "control")

# Subset new_metadata to include only samples with treatment condition "treated"

new_metadata_treated <- subset(new_metadata, treatment.ch1 == "treated")

# Create new_metadata_susceptible containing only susceptible samples
new_metadata_susceptible <- subset(new_metadata, genotype.ch1 == "susceptible")

# Create new_metadata_resistant containing only resistant samples
new_metadata_resistant <- subset(new_metadata, genotype.ch1 == "resistant")





# Step 2: construct a DESeqDataSet object ----------

asd <- DESeqDataSetFromMatrix(countData = raw_data_matrix,
                              colData = new_metadata,
                              design = ~ genotype.ch1 + treatment.ch1
                              )
asd

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(asd)) >= 10
keep
asd <- asd[keep,]
asd

# set the factor level

asd$treatment.ch1 <- relevel(asd$treatment.ch1, ref = "control")



# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
asd <- DESeq(asd)
# Extract differential expression results
resT <- results(asd)
# Display the results
resT
# Summarize the results
summary(resT)
# Write the results to a file
dgh <- write.table(resT, file = "res_deg_treated_vs_control.xls", sep = "\t", quote = FALSE, row.names = TRUE)










#Control Susceptible vs Control Resistant
#Construct a DESeqDataSet object

asdC <- DESeqDataSetFromMatrix(countData = data.matrix(control_raw_data),
                               colData = new_metadata_control,
                               design = ~ genotype.ch1)
asdC
# Run DESeq analysis
asdc <- DESeq(asdC)
# Extract differential expression results
resCS<- results(asdc)
# Display the results
resCS
# Summarize the results
summary(resCS)
# Write the results to a file
gdg <- write.table(resCS, file = "res_deg_control_susceptible_vs_control_resistant.xls", sep = "\t", quote = FALSE, row.names = TRUE)

# Filter for significant DEGs
resCS <- as.data.frame(resCS) %>%
  arrange(padj) %>%
  dplyr::filter(abs(log2FoldChange) >= 1, padj < 0.05)
nrow(resCS)

resCS <- resCS %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not significant"
  ))
summary(resCS)
# Counting up-regulated genes
upregulated_count <- sum(resCS$group == "Up-regulated")

# Counting down-regulated genes
downregulated_count <- sum(resCS$group == "Down-regulated")

# Print the counts
cat("Number of up-regulated genes:", upregulated_count, "\n")
cat("Number of down-regulated genes:", downregulated_count, "\n")


#Control Susceptible vs Treated Susceptible
#Construct a DESeqDataSet object

asdS <-  DESeqDataSetFromMatrix(countData = data.matrix(susceptible_raw_data),
                                colData = new_metadata_susceptible,
                                design = ~ treatment.ch1)
asdS
# Run DESeq analysis
asdS <- DESeq(asdS)
# Extract differential expression results
resCSTS <- results(asdS)
# Display the results
resCSTS
# Summarize the results
summary(resCSTS)
# Write the results to a file
fdf <- write.table(resCSTS, file = "res_deg_control_susceptible_vs_treated_susceptible.xls", sep = "\t", quote = FALSE, row.names = TRUE)

# Filter for significant DEGs
resCSTS <- as.data.frame(resCSTS) %>%
  arrange(padj) %>%
  dplyr::filter(abs(log2FoldChange) >= 1, padj < 0.05)
nrow(resCSTS)

resCSTS <- resCSTS %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not significant"
  ))
summary(resCSTS)


# Counting up-regulated genes
upregulated_count <- sum(resCSTS$group == "Up-regulated")

# Counting down-regulated genes
downregulated_count <- sum(resCSTS$group == "Down-regulated")

# Print the counts
cat("Number of up-regulated genes:", upregulated_count, "\n")
cat("Number of down-regulated genes:", downregulated_count, "\n")



##Control Resistant vs Treated Resistant
#Construct a DESeqDataSet object

asdR <- DESeqDataSetFromMatrix(countData = data.matrix(resistant_raw_data),
                               colData = new_metadata_resistant,
                               design = ~ treatment.ch1)

asdR
# Run DESeq analysis
asdR <- DESeq(asdR)
# Extract differential expression results
resCRTR <- results(asdR)
# Display the results
resCRTR
# Summarize the results
summary(resCRTR)
# Write the results to a file
jdh <- write.table(resCRTR, file = "res_deg_control_resistant_vs_treated_resistant.xls", sep = "\t", quote = FALSE, row.names = TRUE)


# Filter for significant DEGs
resCRTR <- as.data.frame(resCRTR) %>%
  arrange(padj) %>%
  dplyr::filter(abs(log2FoldChange) >= 1, padj < 0.05)
nrow(resCRTR)

resCRTR <- resCRTR %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not significant"
  ))
summary(resCRTR)


# Counting up-regulated genes
upregulated_count <- sum(resCRTR$group == "Up-regulated")

# Counting down-regulated genes
downregulated_count <- sum(resCRTR$group == "Down-regulated")

# Print the counts
cat("Number of up-regulated genes:", upregulated_count, "\n")
cat("Number of down-regulated genes:", downregulated_count, "\n")





#Treated resistant vs treated susceptible
#Construct a DESeqDataSet object

asdT <- DESeqDataSetFromMatrix(countData = data.matrix(treated_raw_data),
                               colData = new_metadata_treated,
                               design = ~ genotype.ch1)
asdT
# Run DESeq analysis
asdT <- DESeq(asdT)
# Extract differential expression results
resTRTS <- results(asdT)
# Display the results
resTRTS
# Summarize the results
summary(resTRTS)
# Write the results to a file
kdd <- write.table(resTRTS, file = "res_deg_treated_resistant_vs_treated_susceptible.xls", sep = "\t", quote = FALSE, row.names = TRUE)


# Filter for significant DEGs
resTRTS <- as.data.frame(resTRTS) %>%
  arrange(padj) %>%
  dplyr::filter(abs(log2FoldChange) >= 1, padj < 0.05)
nrow(resTRTS)

resTRTS <- resTRTS %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not significant"
  ))
summary(resTRTS)


# Counting up-regulated genes
upregulated_count <- sum(resTRTS$group == "Up-regulated")

# Counting down-regulated genes
downregulated_count <- sum(resTRTS$group == "Down-regulated")

# Print the counts
cat("Number of up-regulated genes:", upregulated_count, "\n")
cat("Number of down-regulated genes:", downregulated_count, "\n")

# Step 4: rlog Transformation ----------------------
# Perform rlog transformation
rld <- rlog(asd, blind=FALSE)

# Step 5: PCA Plot ----------------------
# Extract the rlog-transformed values
rlog_matrix <- assay(rld)

# Calculate the PCA

pca <- prcomp(t(rlog_matrix),center = TRUE, 
              scale. = TRUE)
colData(asd)
summary(pca)
str(pca)

# Combine genotype and treatment into a single descriptive factor
colData(asd)$condition <- factor(paste(colData(asd)$genotype.ch1, colData(asd)$treatment.ch1, sep = "_"))

# Check levels of Condition
levels(colData(asd)$condition)

# Assuming pca$x contains PCA results
pca_data <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Condition = colData(asd)$condition
)
#Creating PCA Plot

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point() +
  labs(title = "PCA of rlog-transformed expression data",
       x = paste0("PC1: ", round(pca$sdev[1] / sum(pca$sdev) * 100, 2), "% Variance"),
       y = paste0("PC2: ", round(pca$sdev[2] / sum(pca$sdev) * 100, 2), "% Variance")) +
  theme_minimal()

#variance stabilizing transformation (VST)

vsd <- vst(asd, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(asd, blind = FALSE)
head(assay(rld), 3)

library("dplyr")
library("ggplot2")

asd <- estimateSizeFactors(asd)

df <- bind_rows(
  as_data_frame(log2(counts(asd, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

#Dendogram
# Load necessary libraries
library(DESeq2)
if (!requireNamespace("ggtree", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ggtree",force = TRUE)
}
install.packages("ape",force = TRUE)
library(ggtree)
library(ape)
# Step 1: Extract normalized counts
normalized_counts <- counts(asd, normalized = TRUE)

# Step 2: Filter by significant DEGs (e.g., adjusted p-value < 0.05)
sig_genes <- rownames(resT)[resT$padj < 0.05]
sig_counts <- normalized_counts[sig_genes, ]

# Step 3: Perform hierarchical clustering (log-transform counts)
log_counts <- log1p(sig_counts)
dist_matrix <- dist(t(log_counts))  # Transpose to cluster samples, not genes
hc <- hclust(dist_matrix)

# Step 4: Plot the dendrogram using ggtree
ggtree(as.phylo(hc)) + 
  geom_tiplab() + 
  theme_tree2() + 
  labs(title = "Sample Clustering Based on DEGs")

#Sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDists

#pheatmap

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype.ch1, vsd$treatment.ch1, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

poisd <- PoissonDistance(t(counts(asd)))

#We plot the heatmap in a Figure below.

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( asd$genotype.ch1, asd$treatment.ch1, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# Use the plotPCA function from the DESeq2 package
pcaPlot <- DESeq2::plotPCA(vsd, intgroup=c("genotype.ch1", "treatment.ch1"))

# Print the PCA plot
print(pcaPlot)

#META ANALYSIS

install.packages("metafor")
library(metafor)

# Extract log2 fold change and p-values from DESeq results
logFC1 <- resCS$log2FoldChange
logFC2 <- resT$log2FoldChange
pvalue1 <- resCS$pvalue
pvalue2 <- resT$pvalue
padj1 <- resCS$padj
padj2 <- resT$padj
# Combine the data into a data frame
meta_data1 <- data.frame(logFC1, logFC2, pvalue1, pvalue2)
# Perform meta-analysis
meta_result <- rma.uni(yi = logFC1, sei = sqrt(logFC1^2 / (4 * pvalue1)), data = meta_data)

###Start after running
nrow(asd)
smallestGroupSize <- 4
keep <- rowSums(counts(asd) >= 10) >= smallestGroupSize
asd <- asd[keep,]
nrow(asd)

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)


log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

#variance stabilizing transformation (VST)

vsd <- vst(asd, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(asd, blind = FALSE)
head(assay(rld), 3)

library("dplyr")
library("ggplot2")

asd <- estimateSizeFactors(asd)

df <- bind_rows(
  as_data_frame(log2(counts(asd, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

sampleDists <- dist(t(assay(vsd)))
sampleDists

#pheatmap

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype.ch1, vsd$treatment.ch1, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

poisd <- PoissonDistance(t(counts(asd)))

#We plot the heatmap in a Figure below.

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( asd$genotype.ch1, asd$treatment.ch1, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = c("genotype.ch1", "treatment.ch1"))

pcaData <- plotPCA(vsd, intgroup = c( "genotype.ch1", "treatment.ch1"), returnData = TRUE)
pcaData




# Scale the data
scaled_data <- scale(pcaData[,1:2])

# Perform PCA
pca_result <- prcomp(scaled_data)

# Extract variance explained by each principal component
variance_explained <- summary(pca_result)$importance[2,]

#### PCA plot (51.5)

# Assuming your PCA data is stored in a data frame called pcaData

# Perform PCA
pca <- prcomp(pcaData[,1:2])  # considering only the first two principal components for visualization

# Variance explained by each principal component
variance_explained <- summary(pca)$importance[2,]

# Plot PCA
pca_df <- as.data.frame(pca$x)
pca_df$group <- pcaData$group

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = "PCA Plot",
       x = paste("PC1 (", round(variance_explained[1]*100, 2), "% variance)", sep = ""),
       y = paste("PC2 (", round(variance_explained[2]*100, 2), "% variance)", sep = "")) +
  theme_minimal()



# Extract log2 fold change and p-values from DESeq results
logFC1 <- resCS$log2FoldChange
logFC2 <- resT$log2FoldChange
pvalue1 <- resCS$pvalue
pvalue2 <- resT$pvalue
padj1 <- resCS$padj
padj2 <- resT$padj
rawpval <- list(pvalue1, pvalue2)
FC <- list(logFC1, logFC2)
adjpval <- list(padj1,padj2)
data(adjpval)

studies <- c("resCS", "resT")

DE <- mapply(adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))

DE <- as.data.frame(DE)
colnames(DE) <- paste("DE", studies, sep = ".")



# Assuming DE is a data frame or matrix
colnames(DE) <- paste("DE", studies, sep = ".")


par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="resCS", xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="resT", xlab="Raw p-values")

filtered <- lapply(adjpval, FUN=function(pval) which(is.na(pval)))
rawpval[[1]][filtered[[1]]]=NA
rawpval[[2]][filtered[[2]]]=NA

par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="resCS",xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="resT", xlab="Raw p-values")


fishcomb <- fishercomb(rawpval, BHth = 0.05)
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method", xlab = "Raw p-values (meta-analysis)")

invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)
hist(invnormcomb$rawpval, breaks=100, col="grey",main="Inverse normal method", xlab = "Raw p-values (meta-analysis)")

DEresults <- data.frame(DE, "DE.fishercomb"=ifelse(fishcomb$adjpval<=0.05,1,0), "DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0))
head(DEresults)






# Install necessary packages if you haven't already
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("factoextra", quietly = TRUE)) {
  install.packages("factoextra")
}

# Load the libraries
library(ggplot2)
library(factoextra)

# Load necessary libraries
library(ggplot2)
library(factoextra)

# Example data loading, ensure your data is in `vsd`
# vsd <- read.csv("path/to/your/data.csv", row.names = 1)

# Check the structure of vsd
str(vsd)

# Remove non-numeric columns if any
vsd <- vsd[, sapply(vsd, is.numeric)]

# Remove rows with NA values
vsd <- na.omit(vsd)

# Perform PCA
data.pca <- prcomp(vsd, scale. = TRUE)

# Extract the PCA results
pca_results <- data.frame(Sample = rownames(vsd),
                          PC1 = data.pca$x[, 1],
                          PC2 = data.pca$x[, 2])

# Plot the PCA results using ggplot2
ggplot(pca_results, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(vjust = -1) +
  labs(x = paste0("PC1: ", round(summary(data.pca)$importance[2, 1] * 100, 2), "%"),
       y = paste0("PC2: ", round(summary(data.pca)$importance[2, 2] * 100, 2), "%")) +
  theme_minimal()















