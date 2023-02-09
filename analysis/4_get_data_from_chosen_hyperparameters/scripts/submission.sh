#!/bin/bash 
#SBATCH --partition=dev
#SBATCH -n 15
#SBATCH --mem=100GB
#SBATCH --out=/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/4_get_data_from_chosen_hyperparameters/output/SLURM/output.txt
#SBATCH --error=/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/4_get_data_from_chosen_hyperparameters/output/SLURM/error.txt

/project/shefflab/rivanna_config/bulker_crates/databio/lab/default/Rscript retrieve_data.R