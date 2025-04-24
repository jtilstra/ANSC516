#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 100		            ##number of requested CPUs
#SBATCH --time=48:00:00             ##time requested
#SBATCH --job-name tilstra          ##name of job
#SBATCH -A microbiome               ##name of the queue you are using.

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/JoeTilstra/proj_data_download


pwd

#step 2, Load Modules

echo "loading biocontianers and qiime2"
ml biocontainers qiime2



echo "dada2"
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 17 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 285 \
  --p-trunc-len-r 208 \
  --p-n-threads 0 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza

echo "dada2 end time"
date +"%d %B %Y %H:%M:%S"

echo "dada2 stats visualization"
qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

echo "feature table and featuredata summaries"
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file Proj_metadata.txt

#generate feature table qzv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

echo "generate tree for phylogenetic diversity analyses"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

echo "alpha rarefication plot to pick subsampling depth for core metrics (A&B diversity)"
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 100000 \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization alpha-rarefaction.qzv

echo "Alpha and Beta diversity analysis"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 11298 \
  --m-metadata-file Proj_metadata.txt \
  --output-dir core-metrics-results

echo "Alpha group significance testting"
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization core-metrics-results/observed-features-group-significance.qzv

qiime diversity alpha \
  --i-table table.qza \
  --p-metric simpson \
  --o-alpha-diversity core-metrics-results/simpson_diversity.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/simpson_diversity.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization core-metrics-results/simpson-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

echo "Beta group significance testing"
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Proj_metadata.txt \
  --m-metadata-column Farm_type \
  --o-visualization core-metrics-results/unweighted-unifrac-Farm-type-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Proj_metadata.txt \
  --m-metadata-column Sample_time \
  --o-visualization core-metrics-results/unweighted-unifrac-Sample-timing-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Proj_metadata.txt \
  --m-metadata-column Sample_time \
  --o-visualization core-metrics-results/weighted-unifrac-Sample-timing-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Proj_metadata.txt \
  --m-metadata-column Sample_time \
  --o-visualization core-metrics-results/bray-curtis-Sample-timing-group-significance.qzv \
  --p-pairwise


#echo "taxaonomic analysis"
qiime feature-classifier classify-sklearn \
  --i-classifier silva-classifier-v3v4.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

#Remove mitochondria and chloroplast

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#make visualization of taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Rename files so that the filtered one is table.qza
mv table.qza table-with-mitochondria-and-chloroplast.qza
mv table-no-mitochondria-no-chloroplast.qza table.qza

create taxa barplot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Proj_metadata.txt \
  --o-visualization taxa-bar-plots.qzv

#Differential abundance with ANCOM-BC
echo "ANCOM-BC testing during-season sampling"
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file Proj_metadata.txt \
  --p-where "[Sample_time] = 'during-season'" \
  --o-filtered-table Sample-during-season-table.qza

qiime composition ancombc \
  --i-table Sample-during-season-table.qza \
  --m-metadata-file Proj_metadata.txt \
  --p-formula 'Farm_type' \
  --o-differentials ancombc-farm-type.qza

qiime composition da-barplot \
  --i-data ancombc-farm-type.qza \
  --p-significance-threshold 0.05 \
  --o-visualization da-barplot-farm-type.qzv

