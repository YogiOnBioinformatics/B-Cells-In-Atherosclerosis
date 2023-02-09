
# title: "Get data from chosen hyperparameters"
# author: Yogindra Raghav 
# contact: yogi@email.virginia.edu 

### Important!!! 
# All packages were installed and ran using this R executable: 
# `/project/shefflab/rivanna_config/bulker_crates/databio/lab/default/R`

# Load Packages

rm(list = ls())
library(readxl)
library(tidyverse)
library(CATALYST)
library(FlowSOM)
library(flowCore)
library(igraph)
library(uwot)
library(leiden)
library(matrixStats)
library(ggplot2)
library(glue)
library(gtools)
library(ggpubr)
library(multcomp)
library(diffcyt)

# Read Tabular Data & Metadata

panel_read = read_csv("/home/jve4pt/B-Cells-In-Atherosclerosis/analysis/2_create_pipeline_metadata/output/panelfile_B_cell.csv")
data_folder = "/home/jve4pt/B-Cells-In-Atherosclerosis/analysis/1_CSV_to_FCS/output"
metadata = read_csv("/home/jve4pt/B-Cells-In-Atherosclerosis/analysis/2_create_pipeline_metadata/output/metafile_B_cell.csv")

# Read FCS Files into SingleCellExperiment (SCE) Format & Cluster

sce = prepData(data_folder, panel = panel_read, md = metadata, features = panel_read$fcs_colname)
sce = CATALYST::cluster(sce, features="type", xdim=10, ydim=10, maxK=20, seed = 1234)


# FlowSOM

type_list = type_markers(sce)
sce_flowflame = sce2fcs(sce, split_by = NULL, assay = "exprs")
sce_flowsom = FlowSOM(sce_flowflame,  colsToUse = type_list, nClus = 20, seed = 1234 )


# UMAP Matrix 

matrix_4_umap = t(sce@assays@data$exprs)
matrix_4_umap = matrix_4_umap[,type_markers(sce)]


# UMAP Calculation
set.seed(1234)
uwot_umap_2d = uwot::umap(
  matrix_4_umap,
  n_neighbors = 45,
  n_components = 2,
  min_dist = 0.4,
  n_threads = 15,
  n_sgd_threads=15,
  n_epochs = 500,
  ret_nn = TRUE,
  ret_extra = "fgraph",
  nn_method="fnn"
)

uwot_graph = graph_from_adjacency_matrix(
  uwot_umap_2d$fgraph, 
  weighted = "weight", 
  mode = "directed"
)

# Leiden Clustering 

partition = leiden(
  uwot_graph,
  partition_type = "ModularityVertexPartition",
  n_iterations = -1,
  weights = E(uwot_graph)$weight,
  resolution_parameter = 0.6,
  seed = 1234
)

# Transfer Leiden clusters to FlowSOM meta-clusters

sce@colData$cluster_id = partition

# Output Data 

write.table(
  t(sce@assays@data$exprs), 
  file = "/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/4_get_data_from_chosen_hyperparameters/output/data/arcsinh_transformed_matrix.csv", 
  quote=FALSE, 
  sep=",", 
  row.names=FALSE
)

write.table(
  sce@colData, 
  file="/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/4_get_data_from_chosen_hyperparameters/output/data/metadata.csv", 
  quote=FALSE, 
  sep=",", 
  row.names=FALSE
)
