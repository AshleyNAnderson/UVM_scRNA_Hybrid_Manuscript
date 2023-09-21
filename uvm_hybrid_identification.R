# The background code to provide some complicated functions for the same named Rmd file.
library(Seurat)
library(ggplot2)
# For theme_tufte theme
library(ggthemes)
library(gridExtra)
library(patchwork)
library(scDblFinder)
# To install this, see https://github.com/chris-mcginnis-ucsf/DoubletFinder/tree/master/R
library(DoubletFinder)
# For heat map plot
library(ComplexHeatmap)
# For colors
library(circlize)
library(RColorBrewer)
library(ggrepel)

# To make sure the following python setting work, don't check automatically
# load python in the RStudio's preference. If it doesn't work, try to install
# the latest version of RStudio.
# To use some python package installed for scanpy
# See https://support.rstudio.com/hc/en-us/articles/360023654474-Installing-and-Configuring-Python-with-RStudio
library(reticulate)
# Update the following two paths to your own python env and git repo
Sys.setenv(RETICULATE_PYTHON = "/Users/wug/miniconda3/envs/scanpy/bin/python")
# Make sure the above setting is correct so that we can use scanpy
# Apparently there is no need to call the following if the above works.
# use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
sc <- import('scanpy')

# Working directory is the github repo
setwd("/Users/wug/git/Uveal-Melanoma-Hybrid-Cells")
# Make sure the results are reproducible
random_seed <- 123456
set.seed(random_seed)


# Sample IDs
samples <- c(
    # Class 1 Primary
    "UMM062", "UMM065", 
    # Class 2 Primary
    "UMM059", "UMM061", "UMM063", "UMM064", "UMM066", "UMM069",
    # Class 1 Met
    "BSSR0022",
    # Class 2 Met
    "UMM041L", "UMM067L"
)

# First three are tumors and last three mac
tumor_mac_markers <- c("MLANA", "MITF", "DCT", "CD68", "CD163", "CD14")

# Paths to counts. Please update for your directory structure.
# Data can be donwloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139829
dir_name <- "seq_data"

# Dataframe with provided cell type annotation and doubletFinder output
metadata <- read.csv(file = "scUM_NComms_metadata.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE)

# Add formatted barcodes, tumor cell type labelling
metadata$barcode <- paste0(metadata$orig.ident,
                           "_",
                           stringr::str_extract(metadata$X, ".+(?=-)"))
metadata[metadata$CellType == "Class 1A Primary Tumor Cells", "CellType"] <-
    "Class 1A PRAME- Tumor Cells"
metadata[metadata$CellType == "Class 1B PRAME+ Met Tumor Cells", "CellType"] <-
    "Class 1B PRAME+ Tumor Cells"
metadata[metadata$CellType == "Class 2 PRAME- Primary Tumor Cells","CellType"] <-
    "Class 2 PRAME- Tumor Cells"
metadata[metadata$CellType == "Class 2 PRAME+ Primary Tumor Cells","CellType"] <-
    "Class 2 PRAME+ Tumor Cells"

# Matrix to store final cell counts
cell_counts = matrix(nrow = 7, ncol = 11)
rownames(cell_counts) = c("All Cells", "All Tumor", "All Mac", 
                          "Tumor Hybrid 10-25th", "Tumor Hybrid 25th", 
                          "Mac Hybrid 10-25th", "Mac Hybrid 25th")
colnames(cell_counts) = samples

# Plotlist to store sample plots
plotlist <- list(
    "UMM062" = NA,
    "UMM065" = NA,
    "UMM059" = NA,
    "UMM061" = NA,
    "UMM063" = NA,
    "UMM064" = NA,
    "UMM066" = NA,
    "UMM069" = NA,
    "BSSR0022" = NA,
    "UMM041L" = NA,
    "UMM067L" = NA
)



# Load the seurat object for a specific sample. The function tries to load the 
# cachced object first. If it cannot find a cached object, it will try to load 
# from the original matrix files in the passed dir. There is no any pre-processing
# done in this function exception the min.cellls filtering.
load_seurat_object <- function(sample, 
                               dir_name,
                               min.cells = 3,
                               need_transform = FALSE,
                               which.cells = NA, # If this is NA, all cells should be loaded
                               cache_dir = 'cache') {
    cache_file <- paste(cache_dir, paste("seurat_obj_original_", sample, ".rds", sep = ""), sep = "/")
    if (file.exists(cache_file)) {
        seurat_obj = readRDS(cache_file)
        return(seurat_obj)
    }
    # Setup seurat object  
    dir_name <- paste(dir_name, sample, sep = "/")
    umm <- Read10X(data.dir = dir_name,
                   gene.column = 2, # Use gene names instead of ENSEMBL for easy downstream analysis
                   strip.suffix = TRUE)
    colnames(umm) <- paste(sample, colnames(umm), sep = "_")
    if (!is.na(which.cells)) {
        umm <- umm[, which(colnames(umm) %in% which.cells)]
    }
    umm <- CreateSeuratObject(counts = umm,
                              min.cells = min.cells)
    all_genes <- rownames(umm)
    mt_genes <- all_genes[startsWith(all_genes, 'MT-')]
    # Print out the number of mt genes. It should be 37.
    cat(paste("Total mt_genes: ", length(mt_genes), ".\n", sep = ""))
    if (!is.na(mt_genes)) {
        # Add percent mitochondrial gene count, doubletFinder Output, and cell type annotation
        mt_genes_subset <- mt_genes[which(mt_genes %in% rownames(umm))]
        umm[["percent.mt"]] <- PercentageFeatureSet(umm,
                                                    features = mt_genes_subset)
    }
    if (need_transform) {
        if (!is.na(mt_genes)) {
            # Run SCTransform Normalization
            umm <- SCTransform(umm,
                               vars.to.regress = c("percent.mt"))
        }
        else {
            # Run SCTransform Normalization
            umm <- SCTransform(umm)
        }
    }
    # Cache it
    if (!file.exists(cache_dir)) {
        dir.create(cache_dir, showWarnings = FALSE)
    }
    saveRDS(umm, cache_file)
    return(umm)
}

get_n_pca <- function(umm) {
    pc_var <- (umm@reductions[["pca"]]@stdev / sum(umm@reductions[["pca"]]@stdev)) * 100
    cumu <- cumsum(pc_var)
    nDims <- which(cumu > 95 & pc_var < 5)[1] # A nice choice
    return(nDims)
}

# Plot the ump by loading a saved RDS file from function identify_hybrids
run_ump <- function(sample, dir.name, save = FALSE) {
    rds.file <- paste(dir.name, "/seurat_obj_", sample, ".rds", sep = "")
    if (!file.exists(rds.file)) {
        stop(paste("Cannot find ", rds.file, ". Run identify_hybrids first."))
    }
    umm <- readRDS(rds.file)
    # The code below was copied from Sid's simulate_doublets.R
    umm <- RunPCA(umm, verbose = FALSE)
    nDims <- get_n_pca(umm)
    umm <-  FindNeighbors(umm, dims = 1:nDims)
    # Parameters are based on https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html#cluster-the-cells
    umm <- FindClusters(umm, algorithm = 4, verbose = TRUE) # Use leiden algorithm
    # Specify to use Python's implementation. See: https://rdrr.io/github/satijalab/seurat/man/RunUMAP.html
    # Better to use the R-native implementation, which is much faster and produces similar results
    # Use euclidean so that it can be matched with heatmap
    umm <- RunUMAP(umm, dims = 1:nDims, metric = 'euclidean')
    umm <- RunTSNE(umm, aims = 1:nDims) # Just for comparison and run TSNE here too
    # Default plot by coloring using clusters
    print(DimPlot(umm))
    if (save) {
        pdf(paste(dir.name, "/umap_", sample, ".pdf", sep = ""))
        print(LabelClusters(DimPlot(umm), id = 'ident'))
        dev.off()
    }
    return(umm)
}

# Run paga trajectory analysis. Make sure the python path setting correct at the top.
# Call run_umap first before call this function to avoid re-run PCA.
run_paga <- function(seurat_obj,
                     file = NA) {
    # See section 5: https://theislab.github.io/scanpy-in-R/#the-reticulate-approach
    # Convert seurat_obj into Anndata
    adata_seurat <- sc$AnnData(
        # Need to call as.matrix. Apparent returned is dgCMatrix in latest version of Seurat
        X = t(as.matrix(GetAssayData(seurat_obj))), 
        obs = seurat_obj[[]] # Somehow var cannot work as shown in that web site
    )
    adata_seurat$obsm$update(X_umap = Embeddings(seurat_obj, 'umap'))
    adata_seurat$obsm$update(X_pca = Embeddings(seurat_obj, 'pca'))
    # Just re-run neighbor to avoid complicated migrating
    # Try to get the same set of parameters as in run_umap
    # However, the final output may be different from Seurat!
    sc$pp$neighbors(adata_seurat, 
                    random_state = as.integer(random_seed),
                    n_pcs = get_n_pca(seurat_obj))
    sc$tl$paga(adata_seurat, groups = 'seurat_clusters')
    # Get the coordinates
    # For some reason, cannot set random_seed here because of the following error:
    # TypeError: Cannot cast scalar from dtype('float64') to dtype('int64') according to the rule 'safe'
    sc$pl$paga(adata_seurat, 
               plot = FALSE,
               random_state = as.integer(random_seed))
    reset_pos <- reset_paga_pos(seurat_obj)
    # Have to set by using this method from R to Python:
    # http://cran.nexr.com/web/packages/reticulate/vignettes/introduction.html
    paga <- adata_seurat$uns['paga']
    paga$pos <- reset_pos
    adata_seurat$uns$update(paga = paga)
    # Make sure pos is called
    sc$pl$paga(adata_seurat,
               pos = adata_seurat$uns['paga']['pos']$pos,
               save = file)
    return(adata_seurat)
}

# This function is ported from Python with a R implementation
reset_paga_pos <- function(seurat_obj) {
    # Get the sorted clusters
    clusters <- sort(unique(seurat_obj$seurat_clusters))
    umap_median_matrix <- matrix(nrow = length(clusters), ncol = 2)
    # umap positions
    umap_pos <- Embeddings(seurat_obj, 'umap')
    for (cluster in clusters) {
        cluster_cells <- choose_cells(seurat_obj, 'seurat_clusters', cluster)
        cluster_umap_pos <- umap_pos[cluster_cells, ]
        median_pos <- apply(cluster_umap_pos, 2, median)
        umap_median_matrix[as.integer(cluster), ] <- median_pos
    }
    return(umap_median_matrix)
}


# Plot both tum and mac scores with multiple ways
plot_mac_tumor_scores <- function(seurat_obj,
                                  sample = 'Sample',
                                  reduction_alg = "umap",
                                  need_pdf = FALSE) {
    # Make sure the ident
    Idents(seurat_obj) <- 'celltype_orig'
    # Violin plot
    plot <- VlnPlot(seurat_obj, features = c("Tum_Score1", "Mac_Score1"))
    # UMap
    plot <- plot + 
        DimPlot(seurat_obj, group.by = 'doublet', reduction = reduction_alg) + 
        DimPlot(seurat_obj, group.by = 'celltype_hybrid_high_medium', reduction = reduction_alg) +
        FeaturePlot(seurat_obj, features = c("Tum_Score1", "Mac_Score1"), reduction = reduction_alg, pt.size = 0.1)
    print(plot)
}

# Plot heatmap based on a set of features with hierarchical clustering.
# DoHeatMap() in Seurat cannot perform hierachical clustering and provides
# multiple layer annotation for cells
plot_heat_map <- function(seurat_obj, 
                          use_pca = FALSE,
                          gene_features = NA,
                          cells = NA,
                          cell_features = c("seurat_clusters", "celltype_orig", "celltype_hybrid_high_medium", "New_Cell_Type"),
                          file_name = NA) {
    if (use_pca) {
        df.matrix <- t(seurat_obj@reductions$pca@cell.embeddings)
        gene_features <- rownames(df.matrix)
        cells <- colnames(df.matrix)
    }
    else {
        if (is.na(gene_features)) {
            # Want to choose top 200 features for performance reason
            # Returned results are sorted
            gene_features <- VariableFeatures(seurat_obj)[1:200]
        }
        if (is.na(cells)) {
            # Just for test
            cells <- colnames(seurat_obj)
            # cells <- sample(cells, 2000)
        }
        df <- FetchData(seurat_obj, gene_features, cells)
        # Need another way around
        df <- t(df)
        df.matrix <- as.matrix(df)
    }
    # Add annotation
    cell_annotations <- FetchData(seurat_obj, cell_features, cells)
    ta <- columnAnnotation(cluster = cell_annotations[, 1],
                           celltype_orig = cell_annotations[, 2],
                           celltype_hybrid_high_medium = cell_annotations[, 3],
                           new_cell_type = cell_annotations[, 4],
                           col = list(cluster = create_top_colors(cell_annotations[, 1]),
                                      celltype_orig = create_top_colors(cell_annotations[, 2]),
                                      celltype_hybrid_high_medium = create_top_colors(cell_annotations[, 3]),
                                      new_cell_type = create_top_colors(cell_annotations[, 4])))
    map <- Heatmap(df.matrix,
                   show_column_names = FALSE, # No need to plot cell ids
                   show_row_names = length(gene_features) < 50,
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 2.5),
                   clustering_method_rows = "ward.D2",
                   clustering_method_columns = "ward.D2",
                   top_annotation = ta)
    if (is.na(file_name)) {
        print(map)
    }
    else {
        pdf(file_name, width = 28, height = 10)
        print(map)
        dev.off()
    }
}

create_top_colors <- function(categories) {
    # Create colors
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    topics <- sort(unique(categories))
    top.colors <- sample(col_vector, length(topics))
    names(top.colors) <- topics
    return (top.colors)
}

# Automatically annotate cell types based on markers reported in the original paper
# using AddModuleScore.
annotate_cell_types <- function(seurat_obj, n_bins = 24) {
    # Markers used in the original paper
    cell_type_markers <- list(Tumor = c("MLANA", "MITF", "DCT"),
                              T_cell = c("CD3D", "CD3E", "CD8A"),
                              B_cell = c("CD19", "CD79A", "MS4A1"),
                              Plasma_cell = c("IGHG1", "MZB1", "SDC1", "CD79A"),
                              Monocytes_and_macrophages = c("CD68", "CD163", "CD14"),
                              NK_cell = c("FGFBP2", "FCG3RA", "CX3CR1"),
                              #Retinal_pigment_epithelium = c("RPE65"), // Ignore for the time being
                              #Photoreceptor_cell = c("RCVRN"),
                              Fibroblasts = c("FGF7"),
                              Endothelial_cell = c("PECAM1", "VWF"))
    marker_names <- names(cell_type_markers)
    for (marker_name in marker_names) {
        cat(paste("Working on ", marker_name, "...\n", sep = ""))
        seurat_obj <- AddModuleScore(seurat_obj, 
                                     cell_type_markers[marker_name], 
                                     name = marker_name,
                                     nbin = n_bins)
    }
    # Annotate based on the largest scores for individual cells
    cell_type_names <- names(cell_type_markers)
    score_names <- paste(cell_type_names, "1", sep = "")
    scores <- FetchData(seurat_obj, score_names, colnames(seurat_obj))
    max_indices <- apply(scores, 1, which.max)
    cell_types <- cell_type_names[max_indices]
    names(cell_types) <- colnames(seurat_obj)
    seurat_obj[['New_Cell_Type']] <- cell_types
    return (seurat_obj)
}

choose_cells <- function(seurat_obj, 
                         feature_name, # One single feature only 
                         feature_values) { # This may be a list
    feature_data <- FetchData(seurat_obj, feature_name)
    cells <- colnames(seurat_obj)[feature_data[, 1] %in% feature_values]
    return(cells)
}

create_score_plot <- function(seurat_obj,
                              cells1,
                              cells2,
                              score_name,
                              title,
                              colors) {
    # Get 10th and 25th percentile tumor scores
    scores <- FetchData(seurat_obj, paste(score_name, "1", sep = "")) # A DataFrame has one col only
    # Scores for the first set cells only
    scores_1 <- scores[cells1, ]
    scores_2 <- scores[cells2, ]
    # Want to get the numbers
    score_25th <- quantile(scores_1, probs = c(0.25))
    score_10th <- quantile(scores_1, probs = c(0.1))
    count_0.25 <- sum(scores_2 >= score_25th)
    count_0.10_0.25 <- sum(scores_2 >= score_10th & scores_2 < score_25th)
    # Scores for two cells
    cells_1_2 <- c(cells1, cells2)
    scores_1_2 <- scores[cells_1_2, ]
    # Plot chart limits
    counts_1_2 <- seurat_obj$nCount_RNA[cells_1_2]
    celltypes_1_2 <- seurat_obj$celltype_orig[cells_1_2]
    xMax <- max(counts_1_2)
    # yMin <- min(scores_1_2)
    # yMax <- max(scores_1_2)
    # 
    # Dataframe for plotting
    df <- data.frame(
        cells = cells_1_2,
        score = scores_1_2,
        ncountrna = counts_1_2,
        celltype = celltypes_1_2,
        stringsAsFactors = FALSE
    )
    plot <- ggplot(df, aes(x = ncountrna, y=score, color = celltype)) +
        geom_point(size = 2) +
        scale_color_manual(values = colors) +
        theme_classic() +
        geom_hline(yintercept = quantile(scores_1, c(0.1, 0.25)), 
                   linetype = 'dashed') +
        annotate("text", 
                 x = xMax, 
                 y = quantile(scores_1, c(0.1, 0.25)), 
                 label = as.character(c(paste("10th: ", count_0.10_0.25, sep = ""), 
                                        paste("25th: ", count_0.25, sep = ""))),
                 hjust = 1,
                 vjust = -0.5) +
        xlab("UMI Count") +
        ylab(score_name) + 
        ggtitle(title)
    return(plot)
}

identify_hybrids_2 <- function(seurat_obj,
                               cell_type_1,
                               cell_type_2,
                               cell_type_1_score_name,
                               cell_type_2_score_name,
                               file_name = NA) {
    # Try to catch as many types as possible, e.g. for tumor
    celltype_orig <- unique(seurat_obj$celltype_orig)
    cell_type_1_all <- celltype_orig[grepl(cell_type_1, celltype_orig, ignore.case = TRUE)]
    cell_type_2_all <- celltype_orig[grepl(cell_type_2, celltype_orig, ignore.case = TRUE)]
    
    feature_name <- "celltype_orig"
    cells1 <- choose_cells(seurat_obj, feature_name, cell_type_1_all);
    cells2 <- choose_cells(seurat_obj, feature_name, cell_type_2_all);
    seurat_obj <- add_marker_scores_via_group(seurat_obj,
                                              cells1,
                                              cells2,
                                              cell_type_1_score_name,
                                              cell_type_2_score_name)
    # Need some re-format about the cell type
    celltype_orig <- seurat_obj_scored$celltype_orig
    if (length(cell_type_1_all) > 1) {
        which <- celltype_orig %in% cell_type_1_all
        celltype_orig[which] <- cell_type_1
        seurat_obj_scored[["celltype_orig"]] <- celltype_orig
    }
    if (length(cell_type_2_all) > 1) {
        which <- celltype_orig %in% cell_type_2_all
        celltype_orig[which] <- cell_type_2
        seurat_obj_scored[["celltype_orig"]] <- celltype_orig
    }
    plot <- create_score_plot(seurat_obj_scored,
                              cells1,
                              cells2,
                              cell_type_1_score_name,
                              cell_type_1,
                              c("#438539", "#cfcfcf")); # hard-coded
    plot <- plot + 
        create_score_plot(seurat_obj_scored,
                          cells2,
                          cells1,
                          cell_type_2_score_name,
                          cell_type_2,
                          c("#cfcfcf", "#769cd6"));
    if (is.na(file_name)) {
        print(plot)
    }
    else {
        pdf(file_name, width = 14, height = 7)
        print(plot)
        dev.off()
    }
    return(seurat_obj_scored)
}

add_marker_scores_via_groups <- function(seurat_obj,
                                         cells1,
                                         cells2,
                                         cell_type_1_score_name,
                                         cell_type_2_score_name) {
    markers_1 <- FindMarkers(seurat_obj, 
                             ident.1 = cells1, # Apparently we can pass a set of cells too
                             ident.2 = cells2,
                             only.pos = TRUE);
    # Sort it
    markers_1 <- markers_1[order(markers_1$p_val_adj, 
                                 rev(markers_1$avg_log2FC), 
                                 decreasing = FALSE), ]
    markers_2 <- FindMarkers(seurat_obj,
                             ident.1 = cells2,
                             ident.2 = cells1,
                             only.pos = TRUE)
    markers_2 <- markers_2[order(markers_2$p_val_adj, 
                                 rev(markers_2$avg_log2FC), 
                                 decreasing = FALSE), ]
    # View(markers_1[1:100, ])
    # View(markers_2[1:100, ])
    # Add modules
    seurat_obj_scored <- AddModuleScore(seurat_obj, 
                                        list(rownames(markers_1)[1:50]),
                                        name = cell_type_1_score_name)
    seurat_obj_scored <- AddModuleScore(seurat_obj_scored,
                                        list(rownames(markers_2)[1:50]),
                                        name = cell_type_2_score_name)
    return(seurat_obj_scored)
}

features_scatter_plot <- function(seurat_obj,
                                  feature_1,
                                  feature_2,
                                  cells_1 = NA,
                                  cells_2 = NA,
                                  group_by = 'seurat_clusters') {
    plot <- FeatureScatter(seurat_obj,
                           feature1 = feature_1,
                           feature2 = feature_2,
                           group.by = group_by)
    if (is.na(cells_1) || is.na(cells_2)) {
        print(plot)
    }
    # choose scores
    scores1 <- FetchData(seurat_obj, feature_1)[cells_1, ]
    scores2 <- FetchData(seurat_obj, feature_2)[cells_2, ]
    plot <- plot + geom_hline(yintercept = quantile(scores2, c(0.1, 0.25)), 
                              linetype = 'dashed')
    plot <- plot + geom_vline(xintercept = quantile(scores1, c(0.1, 0.25)),
                              linetype = 'dashed')
    print(plot)
    return(plot)
}

ridge_diff_gene_expression_lot <- function(seurat_obj,
                                           cells1,
                                           cells2,
                                           clusters = NA) {
    markers_1 <- FindMarkers(seurat_obj, 
                             ident.1 = cells1, # Apparently we can pass a set of cells too
                             ident.2 = cells2,
                             only.pos = TRUE);
    # Sort it
    markers_1 <- markers_1[order(markers_1$p_val_adj, 
                                 rev(markers_1$avg_log2FC), 
                                 decreasing = FALSE), ]
    if (is.na(clusters)) {
        idents = NULL
    }
    else {
        idents = clusters
    }
    
    plot1 <- RidgePlot(seurat_obj, 
                       rownames(markers_1)[1:12],
                       ncol = 4, 
                       group.by = 'seurat_clusters',
                       idents = idents)
    print(plot1)
    return(plot1)
}

# Function containing entire patient-level analysis, including plotting and 
# saving seurat objects. Requires declared global variables from previous 
# section and string of sample name as function argument
identify_hybrids <- function(sample, outDir) {
    
    # Subset metadata dataframe to cells from sample
    # This file contains all cells that previously passed QC check, so we will 
    # filter cells in counts table that are present in this file
    meta_sub <- metadata[which(metadata$orig.ident == sample),]
    
    # Setup seurat object  
    umm <- load_seurat_object(sample, 
                              dir_name, 
                              min.cells = 3, 
                              need_transform = TRUE,
                              which.cells = meta_sub$barcode)
    # Add some customized data here
    umm[['doublet']] <-
        setNames(object = meta_sub$DF.classifications_0.25_0.09_4494,
                 nm = meta_sub$barcode)
    umm[['celltype_orig']] <- setNames(object = meta_sub$CellType,
                                       nm = meta_sub$barcode)
    # Generate gene lists by DEG testing, score cells
    # We group all of the four tumor cell types together in this analysis
    umm <- SetIdent(object = umm,
                    value = umm$celltype_orig)
    markers_tum <- FindMarkers(
        umm,
        ident.1 = unique(umm$celltype_orig)[which(grepl("Tumor", unique(umm$celltype_orig)))],
        ident.2 = "Macrophages/Monocytes",
        only.pos = TRUE
    )
    # Need to sort it
    markers_tum <- markers_tum[order(markers_tum$p_val_adj, rev(markers_tum$avg_log2FC), decreasing = FALSE), ]
    markers_mac <- FindMarkers(
        umm,
        ident.1 = "Macrophages/Monocytes",
        ident.2 = unique(umm$celltype_orig)[which(grepl("Tumor", unique(umm$celltype_orig)))],
        only.pos = TRUE
    )
    markers_mac <- markers_mac[order(markers_mac$p_val_adj, rev(markers_mac$avg_log2FC), decreasing = FALSE), ]
    nbin <- 24
    if (sample == 'UMM063') nbin <- 12
    umm <- AddModuleScore(umm, list(rownames(markers_tum)[1:50]), name = "Tum_Score", nbin = nbin)
    umm <- AddModuleScore(umm, list(rownames(markers_mac)[1:50]), name = "Mac_Score", nbin = nbin)
    
    # Get tum/mac cells, and counts
    mac <- WhichCells(object = umm,
                      idents = c("Macrophages/Monocytes"))
    tum <- WhichCells(umm, idents = unique(umm$celltype_orig)[which(grepl("Tumor", unique(umm$celltype_orig)))])
    cell_counts["All Cells", sample] <<- ncol(umm)
    cell_counts["All Mac", sample] <<- length(mac)
    cell_counts["All Tumor", sample] <<- length(tum)
    
    # Subset seruat object
    umm_mac <- subset(umm,
                      cells = mac)
    umm_tum <- subset(umm,
                      cells = tum)
    umm_mac_tum <- subset(umm,
                          cells = c(mac, tum))
    
    # Collapse all tumor cell types into 'Tumor' label 
    temp <- umm_mac_tum$celltype_orig
    temp[which(grepl("Tumor", temp))] <- "Tumor"
    umm_mac_tum[['celltype_collapse']] <- temp
    Idents(umm_mac_tum) <- "celltype_collapse"
    
    # Make named vectors for storing hybrids in object
    temp_hyb_high_only <- umm$celltype_orig
    temp_hyb_high_medium <- umm$celltype_orig
    
    # Get 10th and 25th percentile tumor scores
    scores <- FetchData(umm_tum, "Tum_Score1")$Tum_Score1
    score_25th <- quantile(scores, probs = c(0.25))
    score_10th <- quantile(scores, probs = c(0.1))
    
    # Plot chart limits
    xMax <- max(umm_mac_tum$nCount_RNA)
    yMin <- min(umm_mac_tum$Tum_Score1)
    yMax <- max(umm_mac_tum$Tum_Score1)
    
    # Dataframe for plotting
    df <- data.frame(
        cells = colnames(umm_mac_tum),
        score = umm_mac_tum$Tum_Score1,
        ncountrna = umm_mac_tum$nCount_RNA,
        celltype = umm_mac_tum$celltype_collapse,
        stringsAsFactors = FALSE
    )
    df <- df[rev(order(df$celltype)),]
    
    # Macrophage hybrid plot
    mac_plot <- ggplot(df, aes(x = ncountrna, y=score, color = celltype)) +
        geom_point(size = 2) +
        scale_color_manual(values = c("#438539", "#cfcfcf")) +
        theme_classic() +
        geom_hline(yintercept = quantile(scores, c(0.1, 0.25)), 
                   linetype = 'dashed') +
        annotate("text", 
                 x = xMax, 
                 y = quantile(scores, c(0.1, 0.25)), 
                 label = as.character(c("10th", "25th")),
                 hjust = 1,
                 vjust = -0.5) +
        xlab("UMI Count") +
        ylab("Tumor Score") + 
        ggtitle("Macrophages/Monocytes")
    
    # Get Mac Hybrid Counts
    cell_counts["Mac Hybrid 10-25th", sample] <<- length(which(umm_mac$Tum_Score1 >= score_10th & umm_mac$Tum_Score1 < score_25th))
    cell_counts["Mac Hybrid 25th", sample] <<- length(which(umm_mac$Tum_Score1 >= score_25th))
    temp_hyb_high_only[names(which(umm_mac$Tum_Score1 >= score_25th))] <- "Hybrid"
    temp_hyb_high_medium[names(which(umm_mac$Tum_Score1 >= score_25th))] <- "Hybrid High"
    temp_hyb_high_medium[names(which(umm_mac$Tum_Score1 >= score_10th & umm_mac$Tum_Score1 < score_25th))] <- "Hybrid Medium"
    
    # Get 10th and 25th percentile mac scores
    scores <- FetchData(umm_mac, "Mac_Score1")$Mac_Score1
    score_25th <- quantile(scores, probs = c(0.25))
    score_10th <- quantile(scores, probs = c(0.1))
    
    # Plot chart limits
    xMax <- max(umm_mac_tum$nCount_RNA)
    yMin <- min(umm_mac_tum$Mac_Score1)
    yMax <- max(umm_mac_tum$Mac_Score1)
    
    # Dataframe for plotting
    df <- data.frame(
        cells = colnames(umm_mac_tum),
        score = umm_mac_tum$Mac_Score1,
        ncountrna = umm_mac_tum$nCount_RNA,
        celltype = umm_mac_tum$celltype_collapse,
        stringsAsFactors = FALSE
    )
    df <- df[order(df$celltype),]
    
    # Tumor hybrid plot
    tum_plot <- ggplot(df, aes(x = ncountrna, y=score, color = celltype)) +
        geom_point(size = 2) +
        scale_color_manual(values = c("#cfcfcf", "#769cd6")) +
        theme_classic() +
        geom_hline(yintercept = quantile(scores, c(0.1, 0.25)), 
                   linetype = 'dashed') +
        annotate("text", 
                 x = xMax, 
                 y = quantile(scores, c(0.1, 0.25)), 
                 label = as.character(c("10th", "25th")),
                 hjust = 1,
                 vjust = -0.5) +
        xlab("UMI Count") +
        ylab("Mac Score") + 
        ggtitle("Tumor")
    
    # Get Tum Hybrid Counts
    cell_counts["Tumor Hybrid 10-25th", sample] <<- length(which(umm_tum$Mac_Score1 >= score_10th & umm_tum$Mac_Score1 < score_25th))
    cell_counts["Tumor Hybrid 25th", sample] <<- length(which(umm_tum$Mac_Score1 >= score_25th))
    temp_hyb_high_only[names(which(umm_tum$Mac_Score1 >= score_25th))] <- "Hybrid"
    temp_hyb_high_medium[names(which(umm_tum$Mac_Score1 >= score_25th))] <- "Hybrid High"
    temp_hyb_high_medium[names(which(umm_tum$Mac_Score1 >= score_10th & umm_tum$Mac_Score1 < score_25th))] <- "Hybrid Medium"
    
    
    # Store plots, add hybrids to sample metadata, save Seurat Object
    plotlist[[sample]] <<- mac_plot + tum_plot
    umm[["celltype_hybrid_high_only"]] <-  temp_hyb_high_only
    umm[["celltype_hybrid_high_medium"]] <- temp_hyb_high_medium
    saveRDS(umm, paste0(outDir, "/seurat_obj_", sample, ".rds"))
    return(umm)
}

doublets_by_simulation <- function(seurat_obj){
    # Durante doublet scoring
    # TODO: Need to make sure the selected PC numbers are the same as used in other places (e.g. UMAP)
    seurat_obj <- doubletFinder_v3(seurat_obj, 
                                   PCs = 1:10, pN = 0.25, pK = 0.09, 
                                   nExp = 4494, reuse.pANN = FALSE, 
                                   sct = TRUE)
    # Alternate simulated doublet scoring
    sce <- as.SingleCellExperiment(seurat_obj)
    dens <- computeDoubletDensity(sce)
    # sce$DoubletScore <- dens
    # Directly add the new information instead of converting to avoid
    # data loss during converting
    seurat_obj[['DoubletScore']] <- dens
    return(seurat_obj)
}

# This function is based on https://bioconductor.org/books/3.14/OSCA.advanced/doublet-detection.html
doublets_by_cluster <- function(seurat_obj,
                                dir.name){
    # Doublets using clusters
    # Writes results to working directory
    sce <- as.SingleCellExperiment(seurat_obj)
    dbl.out <- findDoubletClusters(sce, clusters = FetchData(seurat_obj, "seurat_clusters")[, 1], get.all.pairs=TRUE)
    write.csv(data.frame(dbl.out[,0:9]), paste(dir.name, "clustered_doublet_scores.csv", sep = "/"))
    write.csv(data.frame(dbl.out[,10]), paste(dir.name, "clustered_doublet_scores_allpairs.csv", sep = "/"))
}

# Do pathway enrichment analysis using Reactome
analyze_pathways <- function(sample_name, file_name, need_pos = T) {
    # Do some filtering
    markers <- read.table(file_name, sep = "\t", header = T)
    if (need_pos) {
        which <- markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.58 # 1.5 fold expression
        markers <- rownames(markers[which, ])
    }
    else {
        which <- markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.58
        markers <- rownames(markers[which, ])
    }
    require(httr)
    # This url is from the Reatome web site
    url <- 'https://reactome.org/AnalysisService/identifiers/projection?includeDisease=false&interactors=false&order=ASC&page=1&pageSize=10000&pValue=1&resource=TOTAL&sortBy=ENTITIES_FDR&species=9606'
    r <- POST(url, body = markers, content_type('text/plain'))
    results <- content(r, 'parsed')
    # Convert results into a Dataframe
    # Note: This is not a R way. Better to use the approach here: 
    # https://www.r-bloggers.com/2018/10/converting-nested-json-to-a-tidy-data-frame-with-r/
    result_df <- data.frame()
    for (i in 1 : length(results$pathways)) {
        pathway <- results$pathways[[i]]
        row <- c(pathway$stId, pathway$dbId, pathway$name, pathway$entities$total,
                 pathway$entities$found, pathway$entities$ratio, pathway$entities$pValue, 
                 pathway$entities$fdr)
        result_df <- rbind(result_df, row)
    }
    colnames(result_df) <- c('stId', 'dbId', 'name', 'total_entities', 'found_entities', 'ratio', 'pValue', 'fdr')
    # Output
    out_file_name <- paste0(strsplit(file_name, '\\.')[[1]][1], '_pos_', need_pos, '_pathways.tsv')
    write.table(result_df, out_file_name, quote = F, sep = '\t', row.names = F)
}

analyze_pathways_all <- function(dir = 'uvm/deg', name_ends_with='_hyb_cluster_vs_mac_011122.tsv') {
    if (!dir.exists(dir)) {
        stop("No directory existing!")
    }
    files <- list.files(dir)
    files <- files[grepl(name_ends_with, files)]
    for (file in files) {
        # Get the sample name from file
        sample = strsplit(file, '_')[[1]][1]
        analyze_pathways(sample, paste0(dir, "/", file), T)
        analyze_pathways(sample, paste0(dir, "/", file), F)
    }
}

# Perform differentila gene expression analysis between the hybrid cluster 
# and tumor cells or macrophage cells. The annotations of cells are based
# the original paper, except cells in the hybrid cluster. 
# The passed uvm should be a Seurat object returned from run_umap().
list_markers <- function(uvm, name, hybrid_cluster, outDir) {
    if (!dir.exists(outDir)) {
        dir.create(outDir)
    }
    # Get the cells in the hyvrid cluster
    hybrid_cells <- choose_cells(uvm, 'seurat_clusters', hybrid_cluster)
    print(paste0("Total hybrid cluster cells: ", length(hybrid_cells)))
    # Get all tumor types based on the following test
    temp <- unique(uvm$celltype_orig)
    tumor_types <- temp[which(grepl("Tumor", temp))]
    tumor_cells <- choose_cells(uvm, 'celltype_orig', tumor_types)
    print(paste0("Total tumor cells: ", length(tumor_cells)))
    # Remove hybrid cells
    shared <- tumor_cells %in% hybrid_cells
    tumor_cells <- tumor_cells[!shared]
    print(paste0("After removing hybrid cluster cells: ", length(tumor_cells)))
    # Get the macrophage cells
    mac_cells <- choose_cells(uvm, 'celltype_orig', 'Macrophages/Monocytes')
    print(paste0("Total mac cells: ", length(mac_cells)))
    shared <- mac_cells %in% hybrid_cells
    mac_cells <- mac_cells[!shared]
    print(paste0("After removing hybrid cluster cells: ", length(mac_cells)))
    # Will generate 4 gene lists:
    # hybrid vs. mac
    markers <- FindMarkers(uvm, 
                           ident.1 = hybrid_cells,
                           ident.2 = mac_cells,
                           only.pos = FALSE); # Keep positive and negative together for easy handling
    # file_name <- paste0(outDir, "/", name, "_hyb_cluster_vs_mac_011122.tsv")
    file_name <- paste0(outDir, "/", name, "_hyb_cluster_vs_mac_012122.tsv")
    write.table(markers, file_name, quote = F, sep = "\t")
    analyze_pathways(markers, file_name, T)
    analyze_pathways(markers, file_name, F)
    # hybrid vs. tumor
    markers <- FindMarkers(uvm, 
                           ident.1 = hybrid_cells,
                           ident.2 = tumor_cells,
                           only.pos = FALSE); # Keeep positive and negative together for easy handling
    # file_name <- paste0(outDir, "/", name, "_hyb_cluster_vs_tumor_011122.tsv")
    file_name <- paste0(outDir, "/", name, "_hyb_cluster_vs_tumor_012122.tsv")
    write.table(markers, file_name, quote = F, sep = "\t")
    analyze_pathways(markers, file_name, T)
    analyze_pathways(markers, file_name, F)
}

create_df_of_scores <- function(dir_name = 'uvm/deg',
                                name_ends_with = '_vs_mac_011122_pos_TRUE_pathways.tsv',
                                is_for_pathway = T) {
    files <- list.files(dir_name)
    files <- files[grepl(name_ends_with, files)]
    # Get the list of all genes or pathways
    genes <- c()
    for (file in files) {
        markers <- read.table(paste0(dir_name, "/", file), sep = "\t", header = T)
        if (is_for_pathway) {
            rownames(markers) <- markers$name
        }
        genes <- c(genes, rownames(markers))
    }
    genes <- unique(genes)
    genes <- sort(genes)
    print(paste0("Total features: ", length(genes)))
    # Create a gene-based matrix
    # This is a named list
    samples <- sapply(files, FUN = function(f) {strsplit(f, '_')[[1]][1]})
    df <- NULL
    for (name in names(samples)) {
        # name actually is a file name
        markers <- read.table(paste0(dir_name, "/", name), sep = "\t", header = T)
        if (is_for_pathway) {
            rownames(markers) <- markers$name
            # Convert to -log10(FDR)
            fdr_score <- -log10(markers$fdr)
        }
        else {
            fdr_score <- -log10(markers$p_val_adj)
        }
        # Get rid of inf and replace it with max + 1
        temp <- fdr_score
        temp[is.infinite(temp)] <- NA
        max <- max(temp, na.rm = T)
        fdr_score[is.infinite(fdr_score)] <- max + 1;
        # sort based on genes
        fdr_score <- fdr_score[match(genes, rownames(markers))]
        if (is.null(df)) {
            df <- data.frame(fdr_score)
        }
        else {
            df <- cbind(df, fdr_score)
        }
    }
    colnames(df) <- samples
    rownames(df) <- genes
    return(df)
}

create_df_of_binary <- function(dir_name = 'uvm/deg',
                                name_ends_with = '_vs_mac_011122_pos_TRUE_pathways.tsv',
                                is_for_pathway = T) {
    files <- list.files(dir_name)
    files <- files[grepl(name_ends_with, files)]
    # Get the list of all genes
    genes <- c()
    for (file in files) {
        markers <- read.table(paste0(dir_name, "/", file), sep = "\t", header = T)
        if (is_for_pathway) {
            which <- markers$fdr < 0.05
            rownames(markers) <- markers$name
        }
        else {
            which <- markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.58 # 1.5 fold expression
        }
        markers <- rownames(markers[which, ])
        genes <- c(genes, markers)
    }
    genes <- unique(genes)
    genes <- sort(genes)
    print(paste0("Total features: ", length(genes)))
    # Create a gene-based matrix
    # This is a named list
    samples <- sapply(files, FUN = function(f) {strsplit(f, '_')[[1]][1]})
    df <- NULL
    for (name in names(samples)) {
        # name actually is a file name
        markers <- read.table(paste0(dir_name, "/", name), sep = "\t", header = T)
        if (is_for_pathway) {
            which <- markers$fdr < 0.05
            rownames(markers) <- markers$name
        }
        else {
            which <- markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.58 # 1.5 fold expression
        }
        markers <- rownames(markers[which, ])
        which <- as.numeric(genes %in% markers)
        if (is.null(df)) {
            df <- data.frame(which)
        }
        else {
            df <- cbind(df, which)
        }
    }
    colnames(df) <- samples
    rownames(df) <- genes
    return(df)
}

plot_sample_gene_heat_map_based_hyb_markers <- function(dir_name, 
                                                        name_ends_with = '_vs_tumor_011122_pos_TRUE_pathways.tsv',
                                                        is_for_pathway = T,
                                                        use_binary = F) {
    df <- NA
    if (use_binary) {
        df <- create_df_of_binary(dir_name, name_ends_with, is_for_pathway)
    }
    else {
        df <- create_df_of_scores(dir_name, name_ends_with, is_for_pathway)
    }
    tumor_type <- c('C1M', 'C2M', 'C2P', 'C2P', 'C2P', 'C1P', 'C2P')
    top_colors <- c('blue', 'red', 'purple', 'yellow')
    names(top_colors) <- sort(unique(tumor_type))
    ta <- columnAnnotation(tumor_type = tumor_type,
                           col = list(tumor_type = top_colors))
    pdf(paste0(dir_name, "/HeatMap_Hyp_vs_Tumor_Shared_Pathways_012422.pdf"), height = 7)
    df <- na.omit(df)
    heatmap <- Heatmap(as.matrix(df),
                       # name = 'Hit',
                       # col = c('white', 'brown'),
                       clustering_distance_rows = ifelse(use_binary, 'binary','euclidean'), 
                       clustering_distance_columns = ifelse(use_binary, 'binary','euclidean'), 
                       top_annotation = ta,
                       clustering_method_rows = "ward.D2",
                       cluster_rows = T,
                       clustering_method_columns = "ward.D2",
                       row_names_gp = gpar(fontsize = 3),
                       show_row_names = T,
                       column_names_gp = gpar(fontsize = 10))
    print(heatmap)
    dev.off()
    # return (df)
}

# plot_sample_gene_heat_map_based_hyb_markers('uvm/deg',
#                                             name_ends_with = '_hyb_cluster_vs_tumor_011122_pos_TRUE_pathways.tsv',
#                                             is_for_pathway = T,
#                                             use_binary = F)

# Run all analysis and export the results into a PDF file
run_all <- function(sample, outDir) {
    # Result dir based on sample
    if (!file.exists(sample)) {
        dir.create(sample)
    }
    file.name <- paste(sample, "/", sample, "_UMAP_Clusters.pdf", sep = "")
    pdf(file.name)
    uvm <- identify_hybrids(sample, outDir)
    uvm <- run_ump(sample, outDir)
    # # Re-annotate the cell types based on gene markers reported in the paper
    uvm <- annotate_cell_types(uvm)
    # Plot cluster labels 
    plot <- DimPlot(uvm)
    # Label clusters
    plot <- LabelClusters(plot = plot, id = 'ident')
    print(plot)
    dev.off()
    file.name <- paste(sample, "/", sample, "_UMAP_CellTypes.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    # Tother with celltype_orig and celltype_hubrid_high_medium
    plot <- DimPlot(uvm) + DimPlot(uvm, group.by = 'celltype_orig') + DimPlot(uvm, group.by = 'celltype_hybrid_high_medium')
    print(plot)
    dev.off()
    # Tumor and Mac score
    file.name <- paste(sample, "/", sample, "_UMAP_Tum_Mac_Scores.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    plot <- DimPlot(uvm) + FeaturePlot(uvm, c("Tum_Score1", "Mac_Score1"), pt.size = 0.1)
    print(plot)
    dev.off()
    # Violin plot of tumor and mac scores
    file.name <- paste(sample, "/", sample, "_Violin_Tum_Mac_Scores.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    plot <- VlnPlot(uvm, features = c('Tum_Score1', 'Mac_Score1'), group.by = 'seurat_clusters') + 
        VlnPlot(uvm, features = c('Tum_Score1', 'Mac_Score1'), group.by = 'celltype_hybrid_high_medium')
    print(plot)
    dev.off()
    # Plot mac and tumor markers
    file.name <- paste(sample, "/", sample, "_Tumor_Mac_Markers.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    markers <- c("MLANA", "MITF", "DCT", "CD68", "CD163", "CD14")
    plot <- DimPlot(uvm, group.by = 'celltype_hybrid_high_medium') + FeaturePlot(uvm, features = markers, pt.size = 0.1)
    print(plot)
    dev.off()
    # marker scores and new cell type
    file.name <- paste(sample, "/", sample, "_UMAP_New_CellTypes.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    plot <- DimPlot(uvm, group.by = 'celltype_orig') + DimPlot(uvm, group.by = 'New_Cell_Type') + 
        FeaturePlot(uvm, c('Tumor1', 'Monocytes_and_macrophages1'), pt.size = 0.1)
    print(plot)
    dev.off()
    # Two heat maps: the legends may be overlap. Plot them manually!!!
    # Heatmap based on PCA
    file.name <- paste(sample, "/", sample, "_HeatMap_PCA.pdf", sep = "")
    plot_heat_map(uvm, use_pca = TRUE, file_name = file.name)
    # Heatmap based on top 200 genes
    file.name <- paste(sample, "/", sample, "_HeatMap_Genes.pdf", sep = "")
    plot_heat_map(uvm, file_name = file.name)
    # Score distributions
    file.name <- paste(sample, "/", sample, "_Score_Distributions.pdf", sep = "")
    identify_hybrids_2(uvm, 
                       "Tumor",
                       "Macrophages/Monocytes",
                       "Tum_Score",
                       "Mac_Score",
                       file.name)
    # Doublet plot
    file.name <- paste(sample, "/", sample, "_UMAP_Doublets.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    plot <- DimPlot(uvm, group.by = 'celltype_orig') + 
        DimPlot(uvm, group.by = 'celltype_hybrid_high_medium') + 
        DimPlot(uvm, group.by = 'doublet')
    print(plot)
    dev.off()
    # Correlation between doublets and gene counts and genes
    file.name <- paste(sample, "/", sample, "_UMAP_Doublet_Counts.pdf", sep = "")
    pdf(file.name, width = 21, height = 7)
    plot <- DimPlot(uvm, group.by = 'doublet') + 
        FeaturePlot(uvm, c("nCount_RNA", "nFeature_RNA"), pt.size =  0.1)
    print(plot)
    dev.off()
    # Simulation-based doublet
    uvm <- doublets_by_simulation(uvm)
    file.name <- paste(sample, "/", sample, "_UMAP_Doublet_Simulations.pdf", sep = "")
    pdf(file.name, width = 14, height = 7)
    plot <- FeaturePlot(uvm, c("pANN_0.25_0.09_4494", 'DoubletScore'), pt.size = 0.1)
    print(plot)
    dev.off()
    # Cluster-based doublet analysis
    doublets_by_cluster(uvm, sample)
}

# To manually check the clusters. This is based on https://satijalab.org/seurat/articles/visualization_vignette.html
check_clusters <- function(seurat_obj) {
    plot <- DimPlot(seurat_obj, group.by = 'seurat_clusters')
    HoverLocator(plot, 
                 information = FetchData(seurat_obj, 
                                         vars = c('seurat_clusters', 
                                                  'celltype_orig')))
}

# Plot genes from differentilal expression analysis for hybrid vs mac and hybrid vs tumor together in 
# the same bar charts. Only rows in the provided list are used
plot_hybrid_vs_mac_vs_tumor_genes <- function(hybrid_vs_tumor_result_file,
                                              hybrid_vs_mac_result_file,
                                              entry_file,
                                              is_for_upgenes = TRUE,
                                              p_value_labels = NULL,
                                              p_value_breaks = NULL) {
    # Load the data files
    hybrid_vs_tumor = read.delim(hybrid_vs_tumor_result_file, sep = '\t')
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    hybrid_vs_mac = read.delim(hybrid_vs_mac_result_file, sep = '\t')
    print(paste0('Size of hybrid_vs_mac: ', nrow(hybrid_vs_mac)))
    entries <- read.csv(entry_file, header = FALSE)
    print(paste0('Total entiries: ', nrow(entries)))
    print("Doing filtering...")
    which <- rownames(hybrid_vs_tumor) %in% entries$V1
    hybrid_vs_tumor <- hybrid_vs_tumor[which, ]
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    which <- rownames(hybrid_vs_mac) %in% entries$V1
    hybrid_vs_mac <- hybrid_vs_mac[which, ]
    print(paste0('Size of hybrid_vs_mac: ', nrow(hybrid_vs_mac)))
    # Create a new dataframe for plot
    hybrid_vs_mac['analysis'] <- 'hybrid vs mac'
    hybrid_vs_mac['-Log10(p_val_adj)'] <- log10(hybrid_vs_mac['p_val_adj'])
    hybrid_vs_tumor['analysis'] <- 'hybrid vs tumor'
    hybrid_vs_tumor['-Log10(p_val_adj)'] <- -log10(hybrid_vs_tumor['p_val_adj'])
    plot_df <- hybrid_vs_mac
    plot_df <- rbind(plot_df, hybrid_vs_tumor)
    plot_df['Gene'] <- rownames(plot_df)
    plot_df <- plot_df[order(plot_df['-Log10(p_val_adj)'], decreasing = F), ]
    # Make sure the col works like the following
    plot_df['Gene'] <- factor(plot_df[, 'Gene'], levels=plot_df[, 'Gene'])
    # Choose either up-regulated genes or download regulated genes
    if (is_for_upgenes) {
        which <- plot_df[, 'avg_log2FC'] > 0
    }
    else {
        which <- plot_df[, 'avg_log2FC'] < 0
    }
    plot_df <- plot_df[which, ]
    View(plot_df)
    # Plot is based on http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
    g_sig <- ggplot(plot_df, aes(x = `Gene`, y = `-Log10(p_val_adj)`, label = analysis)) + 
        geom_bar(stat = 'identity', aes(x = `Gene`, y = `-Log10(p_val_adj)`, fill = analysis), width = 0.65) + 
        scale_fill_manual(name = "Analysis",
                          labels = c("hybrid vs mac", "hybrid vs tumor"), # Make sure the orders are the same between labels and values
                          values = c("hybrid vs mac"="#00ba38", "hybrid vs tumor"="#f8766d")) +
        theme_tufte() +
        theme(axis.text.y.left = element_text(size = 7)) + 
        coord_flip()
    if (p_value_labels && p_value_breaks) {
        g_sig <- g_sig + scale_y_continuous(breaks = p_value_breaks,
                                            labels = p_value_labels)
    }
    g_fold <- ggplot(plot_df, aes(x = `Gene`, y = `avg_log2FC`, label = analysis)) + 
        geom_bar(stat = 'identity', aes(x = `Gene`, y = `avg_log2FC`, fill = analysis), width = 0.65) + 
        scale_fill_manual(name = "Analysis",
                          labels = c("hybrid vs mac", "hybrid vs tumor"),
                          values = c("hybrid vs mac"="#00ba38", "hybrid vs tumor"="#f8766d")) + 
        theme_tufte() +
        theme(axis.text.y.left = element_text(size = 7)) + 
        coord_flip()
    grid.arrange(g_sig, g_fold, ncol = 2)
    return(plot_df)
}

# Plot genes from differentilal expression analysis for hybrid vs mac or hybrid vs tumor in 
# a scatter plot, which basically is simplified volcano plot.
# Plot is based on https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/
scatter_plot_hybrid_vs_mac_or_tumor_genes <- function(hybrid_vs_tumor_result_file,
                                                      main='Main',
                                                      entry_file,
                                                      highlight_genes_file = NULL) {
    # Load the data files
    hybrid_vs_tumor = read.delim(hybrid_vs_tumor_result_file, sep = '\t')
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    entries <- read.csv(entry_file, header = FALSE)
    print(paste0('Total entiries: ', nrow(entries)))
    print("Doing filtering...")
    which <- rownames(hybrid_vs_tumor) %in% entries$V1
    hybrid_vs_tumor <- hybrid_vs_tumor[which, ]
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    hybrid_vs_tumor['-Log10(p_val_adj)'] <- -log10(hybrid_vs_tumor['p_val_adj'])
    # Create a new dataframe for plot
    plot_df <- hybrid_vs_tumor
    plot_df['Gene'] <- rownames(plot_df)
    plot_df <- plot_df[order(plot_df['-Log10(p_val_adj)'], decreasing = F), ]
    plot_df['Expression'] <- apply(plot_df, 1, FUN = function(x) {
        if (as.numeric(x['p_val_adj']) < 0.05) {
            if (as.numeric(x['avg_log2FC']) > 0.5) {
                return('Up Regulated')
            }
            if (as.numeric(x['avg_log2FC']) < -0.5) {
                return('Down Regulated')
            }
        }
        return('Not Significant')
    })
    # Make sure the col works like the following
    plot_df['Gene'] <- factor(plot_df[, 'Gene'], levels=plot_df[, 'Gene'])
    # Genes for highlight
    hilite_df <- plot_df
    if (!is.null(highlight_genes_file)) {
        hilite_genes <- read.csv(highlight_genes_file, header = F)
        which <- rownames(hilite_df) %in% hilite_genes$V1
        hilite_df <- hilite_df[which, ]
    }
    View(plot_df)
    View(hilite_df)
    # Make sure the colors are the same
    expression_type = unique(plot_df[,'Expression'])
    expression_type_colors = c("dodgerblue3", "gray50", "firebrick3")
    if (length(expression_type) == 2) {
        expression_type_colors = c("dodgerblue3", "firebrick3")
    }
    g_sig <- ggplot(plot_df, aes(x = `avg_log2FC`, y = `-Log10(p_val_adj)`)) +
        geom_point(aes(color=Expression), size=2) +
        scale_color_manual(values = expression_type_colors) +
        guides(colour = guide_legend(override.aes = list(size=1.5))) +
        geom_label_repel(hilite_df, mapping = aes(x = `avg_log2FC`, y = `-Log10(p_val_adj)`, label = `Gene`), size = 3, min.segment.length = 0.5) + # Make sure there is a line
        ggtitle(main) + 
        theme_clean() +
        theme(plot.title = element_text(size=14, face = 'bold'),
              plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
              axis.title.x = element_text(size=12, face = 'bold'),
              axis.title.y = element_text(size=12, face = 'bold'))
    # theme(axis.text.y.left = element_text(size = 7))
    print(g_sig)
    return(plot_df)
}

plot_hybrid_vs_mac_vs_tumor_pathways <- function(hybrid_vs_tumor_result_file,
                                                 hybrid_vs_mac_result_file,
                                                 entry_file) {
    # Load the data files
    hybrid_vs_tumor = read.delim(hybrid_vs_tumor_result_file, sep = '\t')
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    hybrid_vs_mac = read.delim(hybrid_vs_mac_result_file, sep = '\t')
    print(paste0('Size of hybrid_vs_mac: ', nrow(hybrid_vs_mac)))
    entries <- read.csv(entry_file, header = FALSE)
    print(paste0('Total entiries: ', nrow(entries)))
    print("Doing filtering...")
    selected_pathways <- NULL
    if (length(colnames(entries)) > 1) {
        selected_pathways = entries$V2
    }
    else{
        selected_pathways = entries$V1    
    }
    which <- tolower(hybrid_vs_tumor$name) %in% tolower(selected_pathways)
    hybrid_vs_tumor <- hybrid_vs_tumor[which, ]
    print(paste0('Size of hybrid_vs_tum: ', nrow(hybrid_vs_tumor)))
    which <- tolower(hybrid_vs_mac$name) %in% tolower(selected_pathways)
    hybrid_vs_mac <- hybrid_vs_mac[which, ]
    print(paste0('Size of hybrid_vs_mac: ', nrow(hybrid_vs_mac)))
    # Create a new dataframe for plot
    hybrid_vs_mac['analysis'] <- 'hybrid vs mac'
    hybrid_vs_mac['-Log10(fdr)'] <- log10(hybrid_vs_mac['fdr'])
    hybrid_vs_tumor['analysis'] <- 'hybrid vs tumor'
    hybrid_vs_tumor['-Log10(fdr)'] <- -log10(hybrid_vs_tumor['fdr'])
    plot_df <- hybrid_vs_mac
    plot_df <- rbind(plot_df, hybrid_vs_tumor)
    plot_df <- plot_df[order(plot_df['-Log10(fdr)'], decreasing = F), ]
    # Make sure the col works like the following
    order <- unique(plot_df$name)
    neworder <- match(tolower(selected_pathways), tolower(order))
    order <- order[neworder]
    plot_df['name'] <- factor(plot_df$name, levels = order)
    # View(plot_df)
    # Plot is based on http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
    g_sig <- ggplot(plot_df, aes(x = `name`, y = `-Log10(fdr)`, label = analysis)) +
        geom_bar(stat = 'identity', aes(x = `name`, y = `-Log10(fdr)`, fill = analysis), width = 0.65) +
        scale_fill_manual(name = "Analysis",
                          labels = c("hybrid vs mac", "hybrid vs tumor"),
                          values = c("hybrid vs mac"="#00ba38", "hybrid vs tumor"="#f8766d")) +
        scale_y_continuous(breaks = c(-10, 10),
                           labels = c(10, 10)) +
        labs(x = 'Pathway') +
        theme_tufte() +
        theme(axis.text.y.left = element_text(size = 10)) +
        coord_flip()
    print(g_sig)
    return(plot_df)
}

# Code for run one sample
# sample <- "BSSR0022"
# sample <- "UMM066"
sample <- "UMM059"
outDir <- "uvm"

# Plot genes together for two differential expression analyses
hybrid_vs_tumor_result_file <- paste0(outDir, '/deg/', sample, '_hyb_cluster_vs_tumor_011122.tsv')
hybrid_vs_tumor_entry_file <- paste0(outDir, '/deg/hybrid_vs_tumor_genes_061023.csv')
hybrid_vs_mac_result_file <- paste0(outDir, '/deg/', sample, '_hyb_cluster_vs_mac_011122.tsv')
hybrid_vs_mac_entry_file <- paste0(outDir, '/deg/hybrid_vs_mac_genes_061023.csv')
# highlite_genes_file <- paste0(outDir, '/deg/GenesToHighlight_052623.txt')
# entry_file <- paste0(outDir, '/deg/DegGeneList_123022.csv')
# highlite_genes_file <- paste0(outDir, '/deg/GenesToHighlight_052623.txt')
# tmp <- plot_hybrid_vs_mac_vs_tumor_genes(hybrid_vs_tumor_result_file,
#                                          hybrid_vs_mac_result_file,
#                                          entry_file)
# tmp <- plot_hybrid_vs_mac_vs_tumor_genes(hybrid_vs_tumor_result_file,
#                                          hybrid_vs_mac_result_file,
#                                          entry_file,
#                                          is_for_upgenes = FALSE,
#                                          p_value_breaks = c(-40, 0, 40, 80),
#                                          p_value_labels = c(40, 0 , 40, 80))
# tmp <- scatter_plot_hybrid_vs_mac_or_tumor_genes(hybrid_vs_tumor_result_file,
#                                                  main='Hybrid vs Tumor',
#                                                  hybrid_vs_tumor_entry_file,
#                                                  hybrid_vs_tumor_entry_file)
# tmp <- scatter_plot_hybrid_vs_mac_or_tumor_genes(hybrid_vs_mac_result_file,
#                                                  main='Hybrid vs Mac',
#                                                  hybrid_vs_mac_entry_file,
#                                                  hybrid_vs_mac_entry_file)

# Plot pathways together
hybrid_vs_tumor_result_file <- paste0(outDir, '/deg/', sample, '_hyb_cluster_vs_tumor_011122_pos_TRUE_pathways.tsv')
hybrid_vs_mac_result_file <- paste0(outDir, '/deg/', sample, '_hyb_cluster_vs_mac_011122_pos_TRUE_pathways.tsv')
# entry_file <- paste0(outDir, '/deg/SigPathwaysList_123022.csv')
entry_file <- paste0(outDir, '/deg/SigPathwaysList_123022_ordered.csv')
tmp <- plot_hybrid_vs_mac_vs_tumor_pathways(hybrid_vs_tumor_result_file,
                                   hybrid_vs_mac_result_file,
                                   entry_file)

# uvm <- run_ump(sample, outDir, save=TRUE)
# uvm_adata <- run_paga(uvm, file = paste("_", sample, '.pdf', sep=""))
# umm <- identify_hybrids(sample, outDir)
# umm <- run_ump(sample, outDir)
# umm <- annotate_cell_types(umm)
# plot_heat_map(umm)
# plot_heat_map(umm, use_pca = TRUE)
# umm <- doublets_by_simulation(umm)
# # Make sure it is saved into the folder having the sample name
# doublets_by_cluster(umm, sample)
# identify_hybrids(sample, outDir)
# identify_hybrids_2(uvm, "T-Cells", "Macrophages/Monocytes", "T_Score", "Mac_Score")
# identify_hybrids_2(uvm, "Tumor", "Macrophages/Monocytes", "Tum_Score", "Mac_Score")