
# title: "Plotting and Clustering"
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


# Iterate through hyperparameters and then create plots

neighbors_parameters = seq(5,100, by=20)
distance_parameters = seq(0.2, 1, by=.2)
resolution_parameters = seq(0.4, .8, by=.2)


for (neighbor_threshold in neighbors_parameters){

  for (distance_threshold in distance_parameters){

    # UMAP Calculation
    set.seed(1234)
    uwot_umap_2d = uwot::umap(
      matrix_4_umap,
      n_neighbors = neighbor_threshold,
      n_components = 2,
      min_dist = distance_threshold,
      n_threads = 40,
      n_sgd_threads=40,
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

    for (resolution_threshold in resolution_parameters){

      output_folder = "/home/jve4pt/B-Cells-In-Atherosclerosis/analysis/3_plotting_and_clustering/output/plots/"
      output_file = glue("neighbors_{neighbor_threshold}_distance_{distance_threshold}_resolution_{resolution_threshold}")

      # Leiden Clustering 

      partition = leiden(
        uwot_graph,
        partition_type = "ModularityVertexPartition",
        n_iterations = -1,
        weights = E(uwot_graph)$weight,
        resolution_parameter = resolution_threshold,
        seed = 1234
      )

      # Transfer Leiden clusters to FlowSOM meta-clusters

      sce@colData$cluster_id = partition
      sce@metadata$cluster_codes$metaL = sce@metadata$cluster_codes$som100
      sce@metadata$cluster_codes$metaL[(max(partition)+1):100] = max(partition)+1

      # UMAP Plotting

      my_color20 = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#469990', '#dcbeff', '#9A6324',  '#800000',  '#808000',  '#2f4f4f', '#a9a9a9', '#ffd8b1','#000000','#fffac8', '#aaffc3', '#f0f8ff', '#FFDAB9', '#A0522D', '#4B0082', '#F5FFFA', '#000080', '#00FFFF')
      umap_data_4_plot = matrix(0, ncol = 3, nrow = length(partition))
      umap_data_4_plot[,1:2] = uwot_umap_2d$embedding
      umap_data_4_plot[,3] =  matrix(partition)

      colnames(umap_data_4_plot) = c("UMAP1", "UMAP2", "leiden")
      umap_data_4_plot = as.data.frame(umap_data_4_plot)
      umap_data_4_plot[,3] = as.factor(umap_data_4_plot[,3])

      umap_figure = ggplot2::ggplot(umap_data_4_plot, aes(x = UMAP1, y = UMAP2)) +
        ggtitle(glue("Neighbors: {neighbor_threshold} | Distance: {distance_threshold} | Resolution: {resolution_threshold}")) +
        geom_point(aes(color = leiden), alpha = .2, size = 0.75) +   
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5), ncol = 1, title = 'Leiden Clusters')) +
        scale_color_manual(values = my_color20) + 
        theme_bw()

      ggsave(glue("{output_file}_UMAP.png"), path=output_folder, umap_figure, height = 6, width = 8, unit = "in")
      
      # Cluster Expression Heatmap 
      set.seed(1234)

      pdf(glue("{output_folder}{output_file}_heatmap.pdf"))
      plot(plotExprHeatmap(sce, features = "type",
                      by = "cluster_id", k = "metaL", scale = "last", k_pal = my_color20,
                      bars = TRUE, perc = TRUE))
      dev.off()

      # Marker Expression Density Plotx
      pdf(glue("{output_folder}{output_file}_density.pdf"))
      plot(plotClusterExprs(sce, k = "metaL", features = "type"))
      dev.off()

    }

  }

}


