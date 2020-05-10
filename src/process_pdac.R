library(data.table)
library(Seurat)
library(optparse)
library(sctransform)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(parallelDist)

self_pca_distances <- function( normal_cells, self_gene_exp_df ) {
    self_normal_dist_mat <- as.matrix(parDist(as.matrix(self_gene_exp_df, dimnames=list(normal_cells,normal_cells))))
    normal_medoid <- which.min(rowSums(self_normal_dist_mat))
    self_dist_normal_from_normal_medoid <- self_normal_dist_mat[,normal_medoid]

    return( self_dist_normal_from_normal_medoid )
}

annotate_normal_edge_cells <- function( malignant_cells, normal_cells, gene_exp_df, normal_cell_type, malignant_cell_type, edge_cell_fraction=0.1, distance="euclidean", self_gene_exp_df=NULL ) {
    all_cells <- c(normal_cells,malignant_cells)
    subset_gene_exp_df <- gene_exp_df[rownames(gene_exp_df) %in% all_cells,]
    dist_mat <- as.matrix(parDist(as.matrix(subset_gene_exp_df,dimnames=list(all_cells,all_cells)),method=distance))

    self_dist_normal_from_normal_medoid <- self_pca_distances( normal_cells, self_gene_exp_df )
    edge_dist_threshold <- quantile( self_dist_normal_from_normal_medoid, 1 - edge_cell_fraction )
    edge_cell_mask <- self_dist_normal_from_normal_medoid > edge_dist_threshold
    cell_category_labels <- rep( "center", length(normal_cells) )
    cell_category_labels[edge_cell_mask] <- "edge"

        
    combined_normal_cells_anno_dt <- data.table()
    mean_dist_normal_to_malignant <- rowMeans(dist_mat[normal_cells,malignant_cells])

    malignant_only_dist_mat <- dist_mat[malignant_cells,malignant_cells]
    normal_only_dist_mat <- dist_mat[normal_cells,normal_cells]
    malignant_medoid <- which.min(rowSums(malignant_only_dist_mat))
    normal_medoid <- which.min(rowSums(normal_only_dist_mat))
    dist_normal_from_malignant_medoid <- dist_mat[normal_cells,malignant_medoid]

    control_idx <- 1
    normal_cells_anno_dt <- data.table()

    dist_normal_from_normal_medoid <- normal_only_dist_mat[,normal_medoid]

    #closer_to_mal_mask <- dist_normal_from_malignant_medoid < mal_dist_threshold
    #$farther_from_mal_mask <- !closer_to_mal_mask 

    new_dt <- data.table( dist_from_normal_medoid=dist_normal_from_normal_medoid, dist_from_malignant_medoid=dist_normal_from_malignant_medoid, cell.name=normal_cells, cell_category=cell_category_labels )
    new_dt$self_dist_from_normal_medoid <- self_dist_normal_from_normal_medoid

    normal_cells_anno_dt <- rbind( normal_cells_anno_dt, new_dt )

    return( normal_cells_anno_dt )
}

read_gene_exp_mat <- function( matrix_path ) {
    #fread will throw a warning that the number of columns it calculated does not match up
    #with the number of columns in line 2. This is because the header is missing a gene name column for
    #the first column of the file. The three lines after the fread call fixes this issue.
    pdac_normal_dt <- fread( matrix_path ) 
    setnames(pdac_normal_dt,names(pdac_normal_dt)[1],"gene_name")
    pdac_col_names <- names(pdac_normal_dt)

    pdac_normal_mat <- as.matrix(pdac_normal_dt[,!c("gene_name")])
    #The following row and column name assignments are needed for input to Seurat

    colnames(pdac_normal_mat) <- pdac_col_names[2:length(pdac_col_names)]
    rownames(pdac_normal_mat) <- pdac_normal_dt$gene_name

    return (pdac_normal_mat)
    #nonzero_counts <- rowSums(pdac_normal_mat > 0)
    #min_expression <- floor(0.1 * ncol(pdac_normal_mat))
    #pdac_normal_mat <- pdac_normal_mat[nonzero_counts > min_expression,1:ncol(pdac_normal_mat)]
}

create_cell_type_list <- function( anno_dt, normal_cell_types_list, malignant_cell_types_list, cell_types_to_exclude_list ) {
    cell_type_list = c()
    if (normal_cell_types_list[1] == "all" && length(normal_cell_types_list) == 1)  {
        normal_cell_types_list = unique(anno_dt$cluster)
    }

    if (malignant_cell_types_list[1] == "all" && length(malignant_cell_types_list) == 1)  {
        malignant_cell_types_list = unique(anno_dt$cluster)
    }

    cell_type_list = unique( c(normal_cell_types_list, malignant_cell_types_list ) )
    cell_type_list = setdiff( cell_type_list, cell_types_to_exclude_list )

    return (cell_type_list)
}

create_full_seurat_object <- function( gene_exp_matrix, anno_dt=NULL ) {
    if (!is.null(anno_dt)) {
        anno_df <- as.data.frame(anno_dt[,!c("cell.name")])
        rownames(anno_df) <- anno_dt$cell.name
        main_pdac_obj <- CreateSeuratObject( gene_exp_matrix, meta.data=anno_df )
    } else {
        main_pdac_obj <- CreateSeuratObject( gene_exp_matrix )
    }

    return (main_pdac_obj)
}

process_subset_with_seurat <- function( main_pdac_obj, anno_dt, cell_type_list  ) {
    subset_anno_dt <- anno_dt[cluster %in% cell_type_list,]
    cell_subset <- subset_anno_dt$cell.name
    seurat_pdac_obj <- subset( main_pdac_obj, cells = cell_subset )
        
    #seurat_pdac_obj <- NormalizeData( seurat_pdac_obj, normalization.method="LogNormalize", scale.factor = 10000 )
    seurat_pdac_obj <- SCTransform( seurat_pdac_obj )
    seurat_pdac_obj <- RunPCA( seurat_pdac_obj, assay="SCT" )
    subset_seurat_obj <- RunUMAP( subset_seurat_obj, reduction="pca", dims=1:50 )

    #seurat_pdac_obj <- FindVariableFeatures(seurat_pdac_obj, selection.method = "vst", nfeatures = 1000, verbose=TRUE )
    #seurat_pdac_obj <- ScaleData(seurat_pdac_obj)

    return (list(seurat_pdac_obj,subset_anno_dt))
}

process_annotations <- function( annotation_path ) {
    pdac_anno_dt <- fread( annotation_path )
    pdac_anno_dt$sample_type <- 'tumour'
    pdac_anno_dt[grepl('^N',cell.name),]$sample_type <- 'normal'
    pdac_anno_dt$sample <- sapply( strsplit( pdac_anno_dt$cell.name, '_' ), function( x ) { return (unlist(x)[[1]]) } )
    pdac_anno_dt$mat_idx <- seq(1,nrow(pdac_anno_dt))

    return (pdac_anno_dt)
}

add_edge_center_annotation <- function( seurat_obj, anno_dt, malignant_cell_types=NULL, normal_cell_types=NULL, num_pcs=50, perform_control = T, 
num_control_shuffles = 100, pairwise_pca=T, no_tumour_adjacent=F ) {
    edge_cells_list = list()
    center_cells_list = list()
    assay_to_use <- "RNA"
    
    DefaultAssay(seurat_obj) <- assay_to_use

    combined_edge_center_dt <- data.table()
    combined_control_dist_dt <- data.table()
    combined_edge_malignant_dist_dt <- data.table()

    if (pairwise_pca == F) {
        seurat_obj <- ScaleData( seurat_obj )
        seurat_obj <- FindVariableFeatures( seurat_obj, selection.method = "vst", nfeatures = 1000, verbose=TRUE,assay=assay_to_use)
        seurat_obj <- RunPCA( seurat_obj, npcs=50, assay=assay_to_use, reduction.name="all_pca", verbose=F)
    }

    loading_mat_list <- list()
    for (normal_cell_type in normal_cell_types) {
        print(normal_cell_type)
        flush.console()
        edge_cells_list[[normal_cell_type]] = list()
        center_cells_list[[normal_cell_type]] = list()

        if (no_tumour_adjacent == TRUE) {
            normal_dt <- anno_dt[cluster == normal_cell_type & sample_type == "normal",]
        } else {
            normal_dt <- anno_dt[cluster == normal_cell_type,]
        }

        #normal_dt <- anno_dt[cluster == normal_cell_type,]
        normal_cells <- normal_dt$cell.name

        if (length(normal_cells) < 50) {
            print("Too few cells available. Skipping")
            next
        }

        normal_subset_seurat_obj <- subset( seurat_obj, cells=normal_cells )
        normal_subset_seurat_obj <- ScaleData( normal_subset_seurat_obj )
        normal_subset_seurat_obj <- FindVariableFeatures( normal_subset_seurat_obj, selection.method = "vst", nfeatures = 1000, verbose=TRUE,assay=assay_to_use)
        normal_subset_seurat_obj <- RunPCA( normal_subset_seurat_obj, npcs=50, assay=assay_to_use, reduction.name="self_pca", verbose=F)
        #loading_mat_list[[normal_cell_type]] <- Loadings(normal_subset_seurat_obj,reduction = "self_pca")
        normal_pca_df <- Embeddings( normal_subset_seurat_obj, reduction="self_pca" )
        self_dist_from_medoid <- self_pca_distances( normal_cells, normal_pca_df )

        if (perform_control == T) {
            control_dist_matrix <- matrix(0,nrow=length(normal_cells),ncol=num_control_shuffles)
            for (i in 1:num_control_shuffles) {
                if (i%%10 == 0) {
                    print(paste0(i,"/",num_control_shuffles))
                    flush.console()
                }
                shuffled_pca_df <- t(apply(normal_pca_df,1,sample))
                self_control_dist_from_medoid <- self_pca_distances( normal_cells, shuffled_pca_df )
                control_dist_matrix[,i] = self_control_dist_from_medoid 
            }
            colnames(control_dist_matrix) <- paste("control",1:num_control_shuffles,sep="_")
            control_dist_dt <- data.table( control_dist_matrix )
            control_dist_dt$cell.name <- normal_cells
            control_dist_dt$actual <- self_dist_from_medoid
            control_dist_dt$cell_type <- normal_cell_type
            combined_control_dist_dt <- rbind( combined_control_dist_dt, control_dist_dt )
        }

        if (length(normal_cells) == 0) {
            print("Skipping")
            next
        }

        for (malignant_cell_type in malignant_cell_types) {
            if (malignant_cell_type == normal_cell_type) {
                next
            }
            #malignant_dt <- anno_dt[cluster == malignant_cell_type & sample_type == "tumour",]
            malignant_dt <- anno_dt[cluster == malignant_cell_type,]
            malignant_cells <- malignant_dt$cell.name

            if (length(malignant_cells) == 0) {
                print("Skipping")
                next
            }

            cells_to_use <- c(malignant_cells,normal_cells)
            subset_seurat_obj <- subset( seurat_obj, cells=cells_to_use )
            if (pairwise_pca) {
                DefaultAssay(subset_seurat_obj) <- assay_to_use

                subset_seurat_obj <- ScaleData( subset_seurat_obj )
                subset_seurat_obj <- FindVariableFeatures( subset_seurat_obj, selection.method = "vst", nfeatures = 1000, verbose=TRUE,assay=assay_to_use)
                var_features <- VariableFeatures( subset_seurat_obj )
                subset_seurat_obj <- DietSeurat( subset_seurat_obj, counts=F, data=T, scale.data = T, features = var_features )
                subset_seurat_obj <- RunPCA( subset_seurat_obj, npcs=50, assay=assay_to_use, reduction.name="pairwise_pca", verbose=F )
                pca_df <- Embeddings( subset_seurat_obj, reduction="pairwise_pca" )
            } else {
                pca_df <- Embeddings( subset_seurat_obj, reduction="all_pca" )
            }

            edge_center_dt <- annotate_normal_edge_cells( malignant_cells, normal_cells, pca_df, self_gene_exp_df=normal_pca_df )
            edge_center_dt$normal_cell_type <- normal_cell_type
            edge_center_dt$malignant_cell_type <- malignant_cell_type
            edge_center_dist_ratio <- mean(edge_center_dt[cell_category == "edge",dist_from_malignant_medoid])/mean(edge_center_dt[cell_category == "center",dist_from_malignant_medoid])
            edge_malignant_dist_dt <- data.table( normal_cell_type = normal_cell_type, malignant_cell_type = malignant_cell_type, edge_center_dist_ratio=edge_center_dist_ratio, dist_type="actual" )

            if (perform_control == T) {
                temp_dt <- copy(edge_center_dt[,.(cell_category, dist_from_malignant_medoid)])
                control_edge_center_dist_ratios <- vector(length=num_control_shuffles)
                for (shuffle_num in 1:num_control_shuffles) {
                    temp_dt$cell_category <- sample(temp_dt$cell_category)
                    control_edge_center_dist_ratios[shuffle_num] <- mean(temp_dt[cell_category == "edge",dist_from_malignant_medoid])/mean(temp_dt[cell_category == "center",dist_from_malignant_medoid])
                }
                control_malignant_dist_dt <- data.table( normal_cell_type=normal_cell_type, malignant_cell_type=malignant_cell_type,edge_center_dist_ratio=control_edge_center_dist_ratios, dist_type="control" )
                edge_malignant_dist_dt <- rbind(edge_malignant_dist_dt,control_malignant_dist_dt) 
            }
            combined_edge_malignant_dist_dt <- rbind( combined_edge_malignant_dist_dt, edge_malignant_dist_dt )
            combined_edge_center_dt <- rbind( combined_edge_center_dt, edge_center_dt )
        }
    }

    #return(loading_mat_list)
    if (perform_control == T) {
        return( list("edge_center_dt"=combined_edge_center_dt,"control_dist_dt"=combined_control_dist_dt,"edge_malignant_dist_dt"=combined_edge_malignant_dist_dt) ) 
    } else {
        return( list("edge_center_dt"=combined_edge_center_dt,"edge_malignant_dist_dt"=combined_edge_malignant_dist_dt) ) 
    }
}

get_entrez_dt <- function( gene_names ) {
    gene_name_to_entrez_df <- select( org.Hs.eg.db, keys=gene_names, columns=c("ENTREZID"), keytype = "SYMBOL" )
    de_entrez_ids <- gene_name_to_entrez_df$ENTREZID

    old_gene_names_vec <- c()
    name_to_entrez_list <- as.list(org.Hs.egALIAS2EG)
    old_gene_names <- gene_name_to_entrez_df[is.na(gene_name_to_entrez_df$ENTREZID),"SYMBOL"]
    entrez_ids_vec <- c()
    for (old_gene_name in old_gene_names) {
        entrez_ids <-  name_to_entrez_list[[old_gene_name]]
        if (length(entrez_ids) >= 1) {
            old_gene_names_vec <- c( old_gene_names_vec, rep( old_gene_name, length(entrez_ids) ) )
            entrez_ids_vec <- c(entrez_ids_vec,entrez_ids)
        } else {
            print(paste("No alias found for",old_gene_name))
        }
    }
    alias_df <- data.frame( SYMBOL=old_gene_names_vec, corrected_entrez_id=entrez_ids_vec )
    gene_name_to_entrez_df <- plyr::join( gene_name_to_entrez_df, alias_df, by="SYMBOL" )

    assign_entrez_id <- function( row ) {

        if (is.na(row[["ENTREZID"]]) && !is.na(row[["corrected_entrez_id"]])) {
            return(row[["corrected_entrez_id"]])
        } else if (is.na(row[["ENTREZID"]]) && is.na(row[["corrected_entrez_id"]])) {
            return(NA)
        } else {
            return(row[["ENTREZID"]])
        }
    }
    gene_name_to_entrez_df$final_entrez_id <- apply( gene_name_to_entrez_df, 1, assign_entrez_id  )
    gene_name_to_entrez_df <- gene_name_to_entrez_df[!is.na(gene_name_to_entrez_df$final_entrez_id),c("SYMBOL","final_entrez_id")]
    entrez_ids <- gene_name_to_entrez_df$final_entrez_id

    return(data.table(gene_name_to_entrez_df))
}

get_grch37_genes_dt <- function(entrez_ids) {
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL",#"ENSEMBL_MART_FUNCGEN", 
                     host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

    grch37_genes_dt <- getBM( attributes=c("chromosome_name","start_position","end_position","strand","entrezgene_id"),#,"go_id"),
                     filters="entrezgene_id", values=entrez_ids, mart=grch37, verbose=F ) %>% data.table(.)

    grch37_genes_dt$chromosome_name <- paste0("chr",grch37_genes_dt$chromosome_name)
    allowed_chroms <- paste0("chr",c(c(1:22,"X","Y")))
    grch37_genes_dt <- grch37_genes_dt[chromosome_name %in% allowed_chroms,]
 
    strand_info <- grch37_genes_dt$strand
    strand_info[strand_info == -1] <- "-"
    strand_info[strand_info == 1] <- "+"
    grch37_genes_dt$strand <- strand_info

    return(grch37_genes_dt)
}

get_edge_de_genes_dt <- function(seurat_obj,edge_center_dt,tumour_adjacent_cells=c()) {
    do_tumour_adjacent_de <- FALSE
    if (length(tumour_adjacent_cells) > 0)
        do_tumour_adjacent_de <- TRUE
        
    normal_cells <- edge_center_dt$cell.name

    edge_ident_df <- data.frame( edge.center.ident=edge_center_dt$cell_category )
    rownames(edge_ident_df) <- normal_cells
    if (do_tumour_adjacent_de) {
        subset_seurat_obj <- subset( seurat_obj, cells=c(normal_cells,tumour_adjacent_cells) )
        tumour_ident_df <- data.frame( edge.center.ident=rep("tumour_adjacent",length(tumour_adjacent_cells)))
        rownames(tumour_ident_df) <- tumour_adjacent_cells
        edge_ident_df <- rbind( edge_ident_df, tumour_ident_df )
        ident_pairs <- list("edge_vs_center"=c("edge","center"),
        "tumour_adjacent_vs_edge"=c("tumour_adjacent","edge"),"tumour_adjacent_vs_center"=c("tumour_adjacent","center"))
    } else {
        subset_seurat_obj <- subset( seurat_obj, cells=normal_cells )
        ident_pairs <- list("edge_vs_center"=c("edge","center"))
    }

    subset_seurat_obj <- AddMetaData( subset_seurat_obj, edge_ident_df, col.name="edge.center.ident" )
    subset_seurat_obj <- SetIdent( subset_seurat_obj, value="edge.center.ident")

    diff_exp_final_dt = data.table()
    for (de_set in names(ident_pairs)) {
        de_pair <- ident_pairs[[de_set]]
        diff_exp_dt <- FindMarkers( subset_seurat_obj, ident.1=de_pair[1],
                                          ident.2=de_pair[2]) %>% data.table(keep.rownames=T) %>% setnames(.,"rn","gene_name")

        diff_exp_dt$gene_type <- "upregulated"
        diff_exp_dt[avg_logFC < 0,]$gene_type <- "downregulated" 
        diff_exp_dt$de_set <- de_set
        diff_exp_final_dt <- rbind(diff_exp_final_dt,diff_exp_dt)
    }

    return(diff_exp_final_dt)
}

create_resampled_seurat_obj <- function( count_mat, meta_data_dt, feature_counts_to_subsample, feature_counts_reference, feature_to_subsample ) {
    set.seed(1024009200)
    cell_names_reference <- names(feature_counts_reference)
    num_cells_to_subsample <- length(feature_counts_to_subsample)
    target_feature_counts <- sort( sample(feature_counts_reference,size=num_cells_to_subsample), decreasing=TRUE )
    feature_counts_to_subsample <- sort(feature_counts_to_subsample,decreasing=TRUE)

    if (feature_to_subsample == "reads")
        resampled_mat <- count_mat[,cell_names_reference]

    cell_num <- 1
    all_genes <- rownames(count_mat)

    cell_names_to_subsample <- names(feature_counts_to_subsample)
    gene_frequencies <- rowMeans(count_mat[,cell_names_to_subsample])
    num_skipped <- 0
    print(dim(count_mat))
    for (cell_name in cell_names_to_subsample) {
        subsampled_feature_count <- target_feature_counts[cell_num]
        if (feature_to_subsample == "genes") {
            genes_expressed_in_cell <- all_genes[count_mat[,cell_name]>0]
            if (subsampled_feature_count > feature_counts_to_subsample[cell_num]) {
                next
            }
            adjusted_gene_frequencies <- gene_frequencies[genes_expressed_in_cell]/sum(gene_frequencies[genes_expressed_in_cell])
            genes_to_keep <- sample( genes_expressed_in_cell, size=subsampled_feature_count, prob=adjusted_gene_frequencies)
            genes_to_remove <- all_genes[!all_genes %in% genes_to_keep]
            count_mat[!all_genes %in% genes_to_keep,cell_name] <- 0
        } else if (feature_to_subsample == "reads") {
            if (subsampled_feature_count > feature_counts_to_subsample[cell_num]) {
                resampled_mat <- RowMergeSparseMatrices( resampled_mat, count_mat[,cell_name] )
                next
            }
            resampled_cell <- SampleUMI( as.matrix(count_mat[,cell_name]), subsampled_feature_count )
            colnames(resampled_cell) <- cell_name
            rownames(resampled_cell) <- all_genes
            resampled_mat <- RowMergeSparseMatrices( resampled_mat, resampled_cell )
        }
        cell_num <- cell_num + 1
    }
    remaining_cells <- meta_data_dt[!cell.name %in% c(cell_names_to_subsample,cell_names_reference),cell.name]
    remaining_mat <- count_mat[,remaining_cells]
    if (feature_to_subsample == "reads") {
        resampled_mat <- RowMergeSparseMatrices( remaining_mat, resampled_mat )
        cell_name_order <- c(remaining_cells,cell_names_reference,cell_names_to_subsample)
    } else if (feature_to_subsample == "genes")  {
        #resampled_mat <- RowMergeSparseMatrices( remaining_mat, count_mat )
        resampled_mat <- count_mat
        cell_name_order <- meta_data_dt$cell.name 
        #print(dim(resampled_mat))
    }
    
    #4 + ""
    new_meta_data_df <- as.data.frame( meta_data_dt[match(cell.name,cell_name_order),.(cell.name,cluster)] )
    rownames(new_meta_data_df) <- new_meta_data_df$cell.name
    resampled_seurat_obj <- CreateSeuratObject( resampled_mat, meta.data=new_meta_data_df )
    new_meta_data_dt <- data.table( resampled_seurat_obj@meta.data, keep.rownames = T ) %>% setnames(.,"rn","cell.name")
    resampled_seurat_obj <- NormalizeData( resampled_seurat_obj )

    return(resampled_seurat_obj)
}

resample_edge_cells <- function( seurat_obj, pdac_anno_dt, edge_center_dt, normal_cell_types, malignant_cell_types, feature_to_subsample, no_tumour_adjacent=T ) {
    meta_data_dt <- data.table( seurat_obj@meta.data, keep.rownames=T ) %>% setnames(.,"rn","cell.name")
    #edge_info <- add_edge_center_annotation( seurat_obj, pdac_anno_dt, malignant_cell_types=malignant_cell_types, normal_cell_types=normal_cell_types, pairwise_pca = T, perform_control=F, no_tumour_adjacent=no_tumour_adjacent )

    for (malignant_cell_type_ in malignant_cell_types) {
        for (normal_cell_type_ in normal_cell_types) {
            if (no_tumour_adjacent == T) {
                sub_meta_data_dt <- meta_data_dt[cluster == malignant_cell_type_ | 
                                                 (cluster == normal_cell_type_ & sample_type == "normal"),]
            } else {
                sub_meta_data_dt <- meta_data_dt[cluster == malignant_cell_type_ | cluster == normal_cell_type_,]
            }

            pair_edge_center_dt <- edge_center_dt[normal_cell_type == normal_cell_type_ & malignant_cell_type == malignant_cell_type_,]
            pair_edge_center_dt <- merge( pair_edge_center_dt, sub_meta_data_dt[,.(cell.name,nCount_RNA,nFeature_RNA)] )
            cell_popn_list = list()
            for (popn_type in c("edge","center")) {
                read_count_df <- pair_edge_center_dt[cell_category == popn_type,.(cell.name,nCount_RNA,nFeature_RNA)]
                if (feature_to_subsample == "reads") {
                    cell_popn_list[[popn_type]] <- read_count_df$nCount_RNA
                } else {
                    cell_popn_list[[popn_type]] <- read_count_df$nFeature_RNA
                }
                names(cell_popn_list[[popn_type]]) <- read_count_df$cell.name
            }

            count_mat <- seurat_obj[["RNA"]]@counts[,sub_meta_data_dt$cell.name]
            resampled_seurat_obj <- create_resampled_seurat_obj( count_mat, sub_meta_data_dt, cell_popn_list[["edge"]], cell_popn_list[["center"]], feature_to_subsample )

            resampled_obj_meta_data_dt <- data.table( resampled_seurat_obj@meta.data )
            resampled_obj_meta_data_dt <- merge( resampled_obj_meta_data_dt, pair_edge_center_dt, by="cell.name", suffixes=c("_after_resampling","_before_resampling"))

            resampled_edge_info <- add_edge_center_annotation( resampled_seurat_obj, pdac_anno_dt, malignant_cell_types=c(malignant_cell_type_),
                normal_cell_types=c(normal_cell_type_), pairwise_pca = T, num_pcs=num_pcs, perform_control=T, no_tumour_adjacent=no_tumour_adjacent )

            resampled_obj_meta_data_dt <- merge( resampled_obj_meta_data_dt, resampled_edge_info$edge_center_dt, by="cell.name",suffixes=c("_after_resampling","_before_resampling")  )
        }
    }

    return(list("resampled_meta_data_dt"=resampled_obj_meta_data_dt,"edge_info"=resampled_edge_info,"seurat_obj"=resampled_seurat_obj))
}
