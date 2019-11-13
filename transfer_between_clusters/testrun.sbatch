#!/bin/bash
#SBATCH --job-name=transfer_between_clusters
#SBATCH --output=transfer_between_clusters.out
#SBATCH --error=transfer_between_clusters.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1

scriptfolder=/project2/haky/yanyul/GitHub/misc-tools/transfer_between_clusters
bash $scriptfolder/transfer.sh /home/yanyul/bin/meme/lib/meme-5.0.1/perl t.cri.yliang@gardner.cri.uchicago.edu:/home/t.cri.yliang/scratch/test/transfer_between_clusters_test $scriptfolder/testrun.log