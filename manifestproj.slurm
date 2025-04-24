#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 40		            ##number of requested CPUs
#SBATCH --time=06:00:00             ##time requested
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


#step 3, delultiplex
echo "manifest import"

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path  projManifest.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path demux.qza


