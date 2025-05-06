install.packages("pheatmap")

# Load the required libraries
library(ggplot2)
library(pheatmap)

data <- read.table("path_abun_unstrat.tsv",
                   header=TRUE, row.names=1, sep="\t")

# Generate a heatmap
pheatmap(data, scale="row", clustering_distance_rows="correlation", 
         clustering_distance_cols="correlation")

#Generate PCA for path abundance
#read in the path abundance file and transpose it
data <- read.table("path_abun_unstrat.tsv", header=TRUE, row.names=1, sep="\t")
data_t <- t(data)
#read in the metadata file
metadata <- read.table("Proj_metadata.txt", header = TRUE, sep = "\t", row.names = 1)
#make a group to sort the points by
group <- metadata$Sample_time 

# Load the data and perform PCA
pca <- prcomp(data_t, scale. = TRUE)
pca_df <- as.data.frame(pca$x[,1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df$group <- metadata[rownames(pca_df), "Sample_time"] 

#plot the PCA of pathway abundance
PCA_Path<- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(title = "PCA of Pathway Abundance",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.title = element_blank())
ggsave("picrust/PCA_of_Pathway_Abundance.png", height = 4, width = 6)

#create DeSeq plot
library(tidyverse)
library(DESeq2)

#load in the data
counts <- read.table("path_abun_unstrat.tsv", header = TRUE, row.names = 1, sep = "\t",
                     check.names = FALSE)
counts <- as.data.frame(t(counts))

#check if sample names match each other
all(rownames(metadata) %in% rownames(counts))  # should return TRUE if not run the code in the 
head(rownames(metadata))                       # in the next few lines to find problem and filter out
head(rownames(counts))
#filter out the NA
metadata <- metadata[colnames(counts), ]
#check to make sure that the sample names match each other 
all(rownames(metadata) %in% rownames(counts)) #should retrun true

#reorder counts in case they don't match metadata file
counts <- counts[rownames(metadata), ]
#transpose Counts
counts <- t(counts)
#verify that everything matches
all(colnames(counts) == rownames(metadata))
#if not
intersect(colnames(counts), rownames(metadata))
setdiff(colnames(counts), rownames(metadata))  # What's in counts but not metadata?
setdiff(rownames(metadata), colnames(counts))  # What's in metadata but not counts?

#remove NA from counts
counts[is.na(counts)] <- 0

#check for missing values in Sample_time and remove them
sum(is.na(metadata$Sample_time))
metadata <- metadata[!is.na(metadata$Sample_time), ]
#create a subset of counts to match the metadata file
counts <- counts[, rownames(metadata)]

#run DeSeq
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Sample_time)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Sample_time", "during.season", "pre.season"))

#view and filter results
res_df <- as.data.frame(res)
res_df$pathway <- rownames(res_df)

#filter significant results
res_sig <- res_df %>% filter(padj < 0.05)

#create a volcano plot of differential pathway abundance
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Differential Pathway Abundance",
       x = "Log2 Fold Change (during vs. pre season)",
       y = "-log10 Adjusted p-value")
ggsave("picrust/Differntial_Pathway_Abundance.png", width = 6, height = 4)




