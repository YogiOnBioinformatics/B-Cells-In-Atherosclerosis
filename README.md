# B Immune Cells in Atherosclerosis

Work done during 3rd UVA BIMS PhD rotation with Professor Coleen McNamara & Stefan Bekiranov (January - February 2023).

# Chronology
1. 📄 `./analysis/1_CSV_to_FCS/scripts/CSV_to_FCS.R`

Convert CSV files to FCS files which are necessary for `CyTOF` analysis. 

2. 📓 `./analysis/2_create_pipeline_metadata/notebooks/create_metadata.ipynb`

Create feature and sample metadata for downstream analysis. 

3. 📄 `./analysis/3_plotting_and_clustering/scripts/plots_and_clusters.R`

Run UMAP and Leiden clustering algorithms, output plots, and do so iteratively across a 3D hyperparameter sweep. 

4. 📄 `./analysis/3_plotting_and_clustering/scripts/submission.sh`

Submit the above `Rscript` for multithreaded analysis. 

5. 📄 `./analysis/4_get_data_from_chosen_hyperparameters/scripts/retrieve_data.R`

Get UMAP and Leiden clustering output with chosen hyperaparameters based on manual validation of outputs from step 3. 

6. 📄 `./analysis/4_get_data_from_chosen_hyperparameters/scripts/submission.sh`

Submit the above `Rscript` for multithreaded analysis. 

7. 📓 `./analysis/5_differential_analysis/notebooks/differential_analysis.ipynb`

Run differential analysis using nonparametric statistics to understand differential abundance of certain proteins in all cells and in a cell-type specific manner, comparing CASE vs CTRL. 
