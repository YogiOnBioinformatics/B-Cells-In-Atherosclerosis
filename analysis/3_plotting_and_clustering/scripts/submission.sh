#!/bin/bash 
#SBATCH --partition=parallel
#SBATCH -N 2
#SBATCH -n 40
#SBATCH --mem=100GB
#SBATCH --out=/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/3_plotting_and_clustering/output/SLURM/output.txt
#SBATCH --error=/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/3_plotting_and_clustering/output/SLURM/error.txt

/project/shefflab/rivanna_config/bulker_crates/databio/lab/default/Rscript plots_and_clusters.R