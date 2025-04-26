
#Install the packages, IF YOU NEED TO :)
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)

#set working directory

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# These files are already in the ANSC516-repo
##############################################

getwd()
###Set your working directory path/to/ANSC516/ANSC-repo/data/moving-pictures
setwd("C:/Users/jatil/OneDrive/Desktop/ANSC516/ANSC516-repo/data/moving-pictures")

list.files()

#This creates a folder called output if one doesn't already exist
if(!dir.exists("output"))
{
  dir.create("output")
}


#Now the qiime2R method
#this is a more automated way to do the same thing as above
metadata<-read_q2metadata("Proj_metadata.txt")
str(metadata)

row.names(metadata) <- metadata$SampleID
metadata <- metadata[,-1]
row.names(metadata)

#reads in the qza files and stores them 
bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
wUF <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")
uUF <- read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")
jac <- read_qza("core-metrics-results/jaccard_pcoa_results.qza")

#makes a way to sepearte body site by colors (i think this is not needed)
body_colors <- c("Red", "Blue")
my_column <- "Sample_time"

#makes a metadata file that is used to plot the ordination
# %>% symbol is used to break up a line of code into multiple lines for tidyverse
#this makes the dataframe for the PC plot
bc_meta <- bc_PCoA$data$Vectors %>%            #goesinto bc_PCoA and selects vectors in the data list
  select(SampleID, PC1, PC2, PC3) %>%              #this selects only the frist 3 PC and the sampleID
  inner_join(metadata, by = c("SampleID" = "SampleID")) #joins the BC PC plots and adds metadata info

#plotting for bray-curtis
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "Sample_time"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() +
  ggtitle("Bray-Curtis Ordination") +
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".png"), height=3, width=4.5, device="png") # save a PDF 3 inches by 4 inches


#weighted unifrac
wUF_meta <- wUF$data$Vectors %>%            #goesinto bc_PCoA and selects vectors in the data list
  select(SampleID, PC1, PC2, PC3) %>%              #this selects only the frist 3 PC and the sampleID
  inner_join(metadata, by = c("SampleID" = "SampleID")) #joins the BC PC plots and adds metadata info

centroids_wuf <- aggregate(cbind(PC1,PC2)~get(my_column),wUF_meta,mean)
colnames(centroids_wuf)[1] <- "Sample_time"

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  ggtitle("Weighted Unifrac Ordination") +
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*wUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*wUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/wUF-ellipse_", my_column,".png"), height=3, width=4.5, device="png")


#unweighted Unifrac
uUF_meta <- uUF$data$Vectors %>%            #goesinto bc_PCoA and selects vectors in the data list
  select(SampleID, PC1, PC2, PC3) %>%              #this selects only the frist 3 PC and the sampleID
  inner_join(metadata, by = c("SampleID" = "SampleID")) #joins the BC PC plots and adds metadata info

centroids_uuf <- aggregate(cbind(PC1,PC2)~get(my_column),wUF_meta,mean)
colnames(centroids_uuf)[1] <- "Sample_time"

ggplot(uUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  ggtitle("Unweighted Unifrac Ordination") +
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*uUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*uUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/uUF-ellipse_", my_column,".png"), height=3, width=4.5, device="png")


#Jaccard
jac_meta <- jac$data$Vectors %>%            #goesinto bc_PCoA and selects vectors in the data list
  select(SampleID, PC1, PC2, PC3) %>%              #this selects only the frist 3 PC and the sampleID
  inner_join(metadata, by = c("SampleID" = "SampleID")) #joins the BC PC plots and adds metadata info

centroids_jac <- aggregate(cbind(PC1,PC2)~get(my_column),wUF_meta,mean)
colnames(centroids_jac)[1] <- "Sample_time"

ggplot(jac_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  ggtitle("Jaccard Ordination") +
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jac$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/jaccard-ellipse_", my_column,".png"), height=3, width=4.5, device="png")




