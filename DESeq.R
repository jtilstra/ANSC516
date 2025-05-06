
## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#to install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


library(qiime2R)
library(phyloseq)
#library(zoo)
library(tidyverse)
library(DESeq2)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
##############################################
getwd
setwd('../data/moving-pictures/')
list.files()

if(!dir.exists("output/taxa"))
  dir.create("output/taxa")

##Qiime2r method of reading in the taxonomy, metadata and table files individually. 
##Run these lines to troubleshoot if you have trouble on line 55.
#taxonomy<-read_qza("taxonomy.qza")
#head(taxonomy$data)
#taxonomy_table<-parse_taxonomy(taxonomy$data)

#metadata<-read_q2metadata("sample-metadata.tsv")
#str(metadata)

#rare_table <- read_qza("core-metrics-results/rarefied_table.qza")
#feature_table <- rare_table$data

##Qiime2R method of creating a phyloseq object
physeq <- qza_to_phyloseq(
  features="core-metrics-results/rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "Proj_metadata.txt"
)

metadata <- data.frame(sample_data(physeq), check.names = F)


#Clean up metadata, slightly
levels(metadata$Sample_time)
metadata$Sample_time.ord = factor(metadata$Sample_time, c("pre.season", "during.season"))
levels(metadata$Sample_time.ord)



#################################################################
###Differential Abundance with DESeq2
#################################################################



#Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.
#If you need help  with DESeq2 install, see this website
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1

#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = FALSE)

OTU.clean2 <- physeq_otu_table + 1


#Now make the phyloseq object:


OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=TRUE)



#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)

#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.

diagdds = phyloseq_to_deseq2(physeq_deseq, ~ Sample_time)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:

#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

alpha = 0.05


run_deseq2 <- function(my_factor, x, y){
  
  my_contrast <- c(my_factor, x, y)
  res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)
  
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
  
  
  #  ###Volcano Plot
  
  #  with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  #  with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  #  with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  
  #Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.
  
  
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  DESeq_fig = ggplot(sigtab, aes(x=Genus, y = log2FoldChange, color=Phylum)) + 
    geom_point(size=3) + 
    ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
    scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)]) +
    #ylim(0,8) +
    geom_text(color="black", x=length(unique(sigtab$Genus))-3, y=max(sigtab$log2FoldChange)-1, label=my_contrast[2], show_guide = F) +
    geom_text(color="black", x=3, y=min(sigtab$log2FoldChange)+1, label=my_contrast[3], show_guide = F) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  
  ggsave(paste0("output/taxa/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 10, width = 35)
}

#call previously defined function with the factors of intrest
run_deseq2("Sample_time", "during.season" , "pre.season")