#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 100                      ##number of requested CPUs
#SBATCH --time=48:00:00             ##time requested
#SBATCH --job-name tilstra          ##name of job
#SBATCH -A microbiome               ##name of the queue you are using.

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/JoeTilstra/proj_data_download


pwd

#step 2, convert qiime files to formats for picurst
qiime tools export \
--input-path rep-seqs.qza \
--output-path rep-seqs

qiime tools export \
--input-path table.qza \
--output-path table.biom

mv ./rep-seqs/* ..
mv ./table.biom/* ..

#step 3, Load Modules

echo "load the picrust conda environment"
conda activate picrust2_env


picrust2_pipeline.py \
  -s filtered_seqs.fna \
  -i feature-table.biom \
  -o picrust2_out_pipeline \
  -p 100 
