library(data.table)
library(Seurat)
library(optparse)
library(sctransform)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library(parallelDist)
library(robustbase)

self_pca_distances <- function( cells, self_pca_df ) {
    self_normal_dist_mat <- as.matrix(parDist(as.matrix(self_pca_df, dimnames=list(cells,cells))))
    medoid <- which.min(rowSums(self_normal_dist_mat))
    self_dist_from_medoid <- self_normal_dist_mat[,medoid]

    return( list("dist_from_medoid"=self_dist_from_medoid,"medoid"=names(medoid)) )
}

annotate_normal_edge_cells <- function( malignant_cells, normal_cells, gene_exp_df, normal_cell_type, malignant_cell_type, edge_cell_fraction=0.1, distance="euclidean", self_gene_exp_df=NULL ) {
    all_cells <- c(normal_cells,malignant_cells)
    subset_gene_exp_df <- gene_exp_df[rownames(gene_exp_df) %in% all_cells,]
    dist_mat <- as.matrix(parDist(as.matrix(subset_gene_exp_df,dimnames=list(all_cells,all_cells)),method=distance))

    ret <- self_pca_distances( normal_cells, self_gene_exp_df )
    self_dist_normal_from_normal_medoid <- ret$dist_from_medoid
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
    dist_normal_from_all_malignant <- rowMeans(dist_mat[normal_cells,malignant_cells])

    control_idx <- 1
    normal_cells_anno_dt <- data.table()

    dist_normal_from_normal_medoid <- normal_only_dist_mat[,normal_medoid]

    #closer_to_mal_mask <- dist_normal_from_malignant_medoid < mal_dist_threshold
    #$farther_from_mal_mask <- !closer_to_mal_mask 

    new_dt <- data.table( dist_from_normal_medoid=dist_normal_from_normal_medoid,
    dist_from_malignant_medoid=dist_normal_from_malignant_medoid, cell.name=normal_cells,
    cell_category=cell_category_labels,dist_normal_from_all_malignant=dist_normal_from_all_malignant )
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

compute_pca <- function( seurat_obj, cells, num_pcs=50, nfeatures=1000, assay="RNA", pcs_to_use=NULL,
vars_to_regress=c()) {
    subset_seurat_obj <- subset( seurat_obj, cells=cells )
    DefaultAssay(subset_seurat_obj) <- assay
    if (assay == "RNA") { 
        subset_seurat_obj <- FindVariableFeatures( subset_seurat_obj, selection.method = "vst", nfeatures = nfeatures, verbose=F )
    }
    subset_seurat_obj <- ScaleData( subset_seurat_obj )
    #subset_seurat_obj <- SCTransform(subset_seurat_obj,verbose=T)
    subset_seurat_obj <- RunPCA( subset_seurat_obj, npcs=num_pcs, verbose=F)
    pca_df <- Embeddings(subset_seurat_obj)
    if (!is.null(pcs_to_use)) {
        pca_df <- pca_df[,pcs_to_use]# %>% tibble::rownames_to_column("cell.name") %>% dplyr::select(c(cell.name,pcs_to_use)) %>%
        #tibble::column_to_rownames("cell.name")
        if (length(pcs_to_use) == 1) {
            pca_df <- as.matrix(pca_df)
        }
        colnames(pca_df) <- pcs_to_use
    }
    pca_loadings <- Loadings(subset_seurat_obj)

    return(list("pca"=pca_df,"loadings"=pca_loadings))
}

select_edge_center_features <- function( seurat_obj, malignant_cell_types=NULL, normal_cell_types=NULL, num_pcs=50, perform_control = T, 
num_control_shuffles = 100, no_tumour_adjacent=F, num_var_features=1000, assay_to_use="RNA", ident_to_use="cluster", edge_cell_fraction=0.1, sample_info_column=NULL ) {
    edge_cells_list = list()
    center_cells_list = list()
    
    DefaultAssay(seurat_obj) <- assay_to_use
    anno_dt <- data.table(seurat_obj@meta.data,keep.rownames=T) %>% setnames("rn","cell.name")

    combined_edge_center_p_value_df <- data.table()
    combined_skewness_p_value_df <- data.table()

    for (normal_cell_type in normal_cell_types) {
        print(normal_cell_type)
        flush.console()
        edge_cells_list[[normal_cell_type]] = list()
        center_cells_list[[normal_cell_type]] = list()

        if (no_tumour_adjacent == TRUE) {
            normal_dt <- anno_dt[get(ident_to_use) == normal_cell_type & sample_type == "normal",]
        } else {
            normal_dt <- anno_dt[get(ident_to_use) == normal_cell_type,]
        }
        normal_cells <- normal_dt$cell.name

        if (length(normal_cells) < 50) {
            print("Too few cells available. Skipping")
            next
        }

        normal_pca <- compute_pca( seurat_obj, cells=normal_cells, num_pcs=num_pcs, nfeatures=num_var_features, assay=assay_to_use)
        ret <- self_pca_distances( cells, normal_pca$pca )
        normal_dist_from_medoid <- ret$dist_from_medoid
        subset_seurat_obj <- subset( seurat_obj, cells=normal_cells )

        #normal_medoid <- ret$medoid
        sample_wise_cell_name_list <- NULL
        if (!is.null(sample_info_column)) {
             samples <- unique(subset_seurat_obj@meta.data[[sample_info_column]]) 
             for (sample_ in samples) { 
                cells <- rownames( subset_seurat_obj@meta.data %>% dplyr::filter(!!sym(sample_info_column) == sample_))
                if (length(cells) == 0)
                    next

                sample_wise_cell_name_list[[sample_]] <- cells
             }
        }

        if (perform_control == T) {
            if (assay_to_use == "RNA") {
                var_features <- subset_seurat_obj %>% NormalizeData %>% FindVariableFeatures(nfeatures=num_var_features) %>% VariableFeatures
                skewness_p_value_df <- perform_gene_shuffle_parallel(normal_cells,subset_seurat_obj[[assay_to_use]]@counts[var_features,],
            num_control_shuffles,sample_wise_cell_name_list,assay_to_use=assay_to_use,npcs=num_pcs) #%>% mutate(normal_cell_type=normal_cell_type)
            } else {
                var_features <- rownames(subset_seurat_obj)
                skewness_p_value_df <- perform_gene_shuffle_parallel(normal_cells,subset_seurat_obj[[assay_to_use]]@scale.data[var_features,],
            num_control_shuffles,sample_wise_cell_name_list,assay_to_use=assay_to_use,npcs=num_pcs) %>% mutate(normal_cell_type=normal_cell_type)
            }
                
            combined_skewness_p_value_df <- rbind(combined_skewness_p_value_df,skewness_p_value_df)
        }

        if (length(normal_cells) == 0) {
            print("Skipping")
            next
        }

        for (malignant_cell_type in malignant_cell_types) {
            if (malignant_cell_type == normal_cell_type) {
                next
            }
            malignant_dt <- anno_dt[get(ident_to_use) == malignant_cell_type,]
            malignant_cells <- malignant_dt$cell.name

            if (length(malignant_cells) == 0) {
                print("Skipping")
                next
            }

            normal_and_malignant_cells <- c(malignant_cells,normal_cells)
            combined_pca <- compute_pca( seurat_obj, cells=normal_and_malignant_cells, num_pcs=num_pcs,nfeatures=num_var_features, assay=assay_to_use )

            for (normal_idx in 1:num_pcs) {
                normal_single_pc_mat <- as.matrix(normal_pca$pca[,normal_idx])
                rownames(normal_single_pc_mat) <- rownames(normal_pca$pca)
                for (combined_idx in 1:num_pcs) {
                    combined_single_pc_mat <- as.matrix(combined_pca$pca[,combined_idx])
                    rownames(combined_single_pc_mat) <- rownames(combined_pca$pca)
                    edge_center_dt <- annotate_normal_edge_cells( malignant_cells, normal_cells, combined_single_pc_mat, normal_cell_type, malignant_cell_type, self_gene_exp_df=normal_single_pc_mat, edge_cell_fraction=edge_cell_fraction )
                    edge_center_dt[,pc:=paste("combined","PC",combined_idx,sep="_")]
                    edge_center_dt$normal_cell_type <- normal_cell_type
                    edge_center_dt$malignant_cell_type <- malignant_cell_type

                    if (perform_control == T) {
                        edge_malignant_dist_dt <- perform_edge_center_shuffle(edge_center_dt,normal_cell_type,malignant_cell_type,num_control_shuffles)
                        edge_malignant_dist_dt[,pc:=paste("combined","PC",combined_idx,sep="_")]
                    }

                    mean_sd_edge_center_df <- edge_malignant_dist_dt %>% dplyr::filter(dist_type != "actual") %>%
                    group_by(normal_cell_type,malignant_cell_type,pc) %>% summarize(m=mean(edge_center_dist_ratio),s=sd(edge_center_dist_ratio))
                    edge_center_p_value_df <- merge( edge_malignant_dist_dt %>% dplyr::filter(dist_type == "actual"),
                                                 mean_sd_edge_center_df, by=c("normal_cell_type","malignant_cell_type","pc")) %>%
                      group_by(normal_cell_type,malignant_cell_type,pc) %>%
                      mutate(p_value=pnorm(edge_center_dist_ratio,m,s),normal_pc=paste0("normal_PC_",normal_idx)) %>% data.table

                    combined_edge_center_p_value_df <- rbind(combined_edge_center_p_value_df,edge_center_p_value_df)
                }
            }
        }
    }

    combined_edge_center_p_value_df <- combined_edge_center_p_value_df %>% mutate(q_value = p.adjust(p_value)) %>% data.table
    return(list("edge_center_p_value"=combined_edge_center_p_value_df,"skewness_p_value"=combined_skewness_p_value_df))
}

perform_gene_shuffle_parallel <- function(normal_cells,gene_exp_df,num_control_shuffles, sample_wise_cell_name_list=NULL, npcs=50,nfeatures=1000, numCores=6,assay_to_use="RNA") {
    shuffled_skewness_mat <- matrix(0,nrow=npcs,ncol=(num_control_shuffles+1))
    features <- rownames(gene_exp_df)
    orig_seurat_obj <- CreateSeuratObject(gene_exp_df) 
    if (assay_to_use == "RNA") { 
        orig_seurat_obj <- NormalizeData(orig_seurat_obj) %>% ScaleData
        scaled_mat <- orig_seurat_obj[["RNA"]]@scale.data
    } else if (assay_to_use == "integrated") {
        scaled_mat <- gene_exp_df
        orig_seurat_obj[["RNA"]]@scale.data <- scaled_mat
    }
    shuffled_seurat_obj <- copy(orig_seurat_obj)

    parallel_pca <- function(idx) {
        shuffled_skewness_vals <- vector(mode="numeric",length=npcs)
        if (idx == 0) {
            obj <- orig_seurat_obj
        } else {
            set.seed(idx)
            if (is.null(sample_wise_cell_name_list)) {
                shuffled_mat <- apply(scaled_mat,1,sample) %>% t
                colnames(shuffled_mat) <- colnames(scaled_mat)
            } else {
                shuffled_mat <- scaled_mat
                for (sample_ in names(sample_wise_cell_name_list)) { 
                    cells <- sample_wise_cell_name_list[[sample_]]
                    if (length(cells) > 1) {
                        shuffled_mat[,cells] <- apply(scaled_mat[,cells],1,sample) %>% t
                    } else {
                        shuffled_mat[,cells] <- scaled_mat[,cells]
                    }
                }
            }
            shuffled_seurat_obj[["RNA"]]@scale.data <- shuffled_mat
            obj <- shuffled_seurat_obj
        }
        print(idx)
        flush.console()
        
        obj <- RunPCA(obj,npcs=npcs,features=features,verbose=F) 
        pca_df <- Embeddings(obj,reduction="pca")
        pc_names <- colnames(pca_df)
        for (pc_idx in seq_along(colnames(pca_df)))  {
            ret <- self_pca_distances( normal_cells, pca_df[,pc_names[pc_idx]] )
            shuffled_skewness_vals[pc_idx] = mc(ret$dist_from_medoid)
        }

        return(shuffled_skewness_vals)
    }

    #skewness_list <- mclapply(0:num_control_shuffles,parallel_pca,mc.cores=numCores)
    skewness_list <- lapply(0:num_control_shuffles,parallel_pca)
    shuffled_skewness_mat <- sapply(skewness_list,as.matrix)

    pc_names <- paste("normal_PC",1:npcs,sep="_")
    rownames(shuffled_skewness_mat) <- pc_names
    colnames(shuffled_skewness_mat) <- c("actual",paste("control",1:num_control_shuffles,sep="_"))

    skewness_df <- reshape2::melt(shuffled_skewness_mat) %>% mutate(Var2=gsub("_[0-9][0-9]*","",Var2))
    mean_sd_skewness_df <- skewness_df %>% dplyr::filter(Var2 != "actual") %>% 
    group_by(Var1) %>% summarize(m=mean(value),s=sd(value))
    skewness_p_value_df <- merge( skewness_df %>% dplyr::filter(Var2 == "actual"),
                                mean_sd_skewness_df, by="Var1") %>% 
     group_by(Var1) %>% mutate(p_value=1-pnorm(value,m,s)) %>% data.frame %>% mutate(q_value=p.adjust(p_value))
    return(skewness_p_value_df)
}

merge_pca_loadings <- function(normal_pca_loadings_mat,combined_pca_loadings_mat) { 
    colnames(combined_pca_loadings_mat) <- paste("combined",colnames(combined_pca_loadings_mat),sep="_")
    colnames(normal_pca_loadings_mat) <- paste("normal",colnames(normal_pca_loadings_mat),sep="_")
    normal_pca_loadings_dt <- data.table( normal_pca_loadings_mat, keep.rownames=T ) %>% setnames("rn","gene")
    combined_pca_loadings_dt <- data.table( combined_pca_loadings_mat, keep.rownames=T ) %>% setnames("rn","gene")
    loadings_dt <- merge(normal_pca_loadings_dt,combined_pca_loadings_dt,by="gene")

    return(loadings_dt)
}

add_edge_center_annotation_classic <- function( seurat_obj, malignant_cell_types=NULL, normal_cell_types=NULL,
num_pcs=50, perform_control = T, num_control_shuffles = 100, no_tumour_adjacent=F, num_var_features=1000,
assay_to_use="RNA", ident_to_use="cluster", edge_cell_fraction=0.1, pooled_pcs_to_use=NULL, normal_pcs_to_use=NULL, sample_info_column=NULL  ) {
    edge_cells_list = list()
    center_cells_list = list()
#    assay_to_use <- "RNA"
    anno_dt <- data.table(seurat_obj@meta.data,keep.rownames=T) %>% setnames("rn","cell.name")
    
    DefaultAssay(seurat_obj) <- assay_to_use

    combined_edge_center_dt <- data.table()
    combined_control_dist_dt <- data.table()
    combined_edge_malignant_dist_dt <- data.table()
    
    if (is.null(normal_pcs_to_use)) {
        normal_pcs_to_use <- paste("normal_PC",1:num_pcs,sep="_")
    }

    if (is.null(pooled_pcs_to_use)) {
        pooled_pcs_to_use <- paste("combined_PC",1:num_pcs,sep="_")
    }

    for (normal_cell_type in normal_cell_types) {
        print(normal_cell_type)
        flush.console()
        edge_cells_list[[normal_cell_type]] = list()
        center_cells_list[[normal_cell_type]] = list()

        if (no_tumour_adjacent == TRUE) {
            normal_dt <- anno_dt[get(ident_to_use) == normal_cell_type & sample_type == "normal",]
        } else {
            normal_dt <- anno_dt[get(ident_to_use) == normal_cell_type,]
        }

        normal_cells <- normal_dt$cell.name

        sample_wise_cell_name_list <- NULL
        if (!is.null(sample_info_column)) {
             samples <- unique(seurat_obj@meta.data[normal_cells,][[sample_info_column]]) 
             for (sample_ in samples) { 
                cells <- intersect( rownames( seurat_obj@meta.data %>% dplyr::filter(!!sym(sample_info_column) == sample_)), normal_cells ) 
                if (length(cells) == 0)
                    next

                sample_wise_cell_name_list[[sample_]] <- cells
             }
        }

        if (length(normal_cells) < 50) {
            print("Too few cells available. Skipping")
            next
        }

        normal_pcs_to_use <- gsub("normal_","",normal_pcs_to_use)
        normal_pca <- compute_pca( seurat_obj, cells=normal_cells, num_pcs=num_pcs, nfeatures=num_var_features, assay=assay_to_use )
        ret <- self_pca_distances( cells, normal_pca$pca[,normal_pcs_to_use] )
        normal_dist_from_medoid <- ret$dist_from_medoid
        if (perform_control == T) {
            if (length(normal_pcs_to_use) > 1 || is.null(normal_pcs_to_use)) {
                control_dist_dt <- perform_pca_shuffle(normal_cells,normal_pca$pca,num_control_shuffles,sample_wise_cell_name_list)
                control_dist_dt[,`:=`(cell.name=normal_cells,actual=normal_dist_from_medoid,cell_type=normal_cell_type)]
                combined_control_dist_dt <- rbind( combined_control_dist_dt, control_dist_dt )
            } else {
                print("Cannot shuffle PCs for skewness test since only a single PC is being used.")
            }
        }

        if (length(normal_cells) == 0) {
            print("Skipping")
            next
        }

        for (malignant_cell_type in malignant_cell_types) {
            if (malignant_cell_type == normal_cell_type) {
                next
            }
            malignant_dt <- anno_dt[get(ident_to_use) == malignant_cell_type,]
            malignant_cells <- malignant_dt$cell.name

            if (length(malignant_cells) == 0) {
                print("Skipping")
                flush.console()
                next
            }

            normal_and_malignant_cells <- c(malignant_cells,normal_cells)
            pooled_pcs_to_use <- gsub("combined_","",pooled_pcs_to_use)
            combined_pca <- compute_pca( seurat_obj, cells=normal_and_malignant_cells, num_pcs=num_pcs,nfeatures=num_var_features, assay=assay_to_use )
            edge_center_dt <- annotate_normal_edge_cells( malignant_cells, normal_cells, as.matrix(combined_pca$pca[,pooled_pcs_to_use]), normal_cell_type, malignant_cell_type, self_gene_exp_df=as.matrix(normal_pca$pca[,normal_pcs_to_use]), edge_cell_fraction=edge_cell_fraction )
            edge_center_dt$normal_cell_type <- normal_cell_type
            edge_center_dt$malignant_cell_type <- malignant_cell_type

            if (perform_control == T) {
                edge_malignant_dist_dt <- perform_edge_center_shuffle(edge_center_dt,normal_cell_type,malignant_cell_type,num_control_shuffles)
                combined_edge_malignant_dist_dt <- rbind( combined_edge_malignant_dist_dt, edge_malignant_dist_dt )
            }
            combined_edge_center_dt <- rbind( combined_edge_center_dt, edge_center_dt )
        }
    }

    merged_pca_loadings_dt <- merge_pca_loadings( normal_pca$loadings, combined_pca$loadings )
    colnames(combined_pca$pca) <- paste("combined",colnames(combined_pca$pca),sep="_")
    colnames(normal_pca$pca) <- paste("normal",colnames(normal_pca$pca),sep="_")
    if (perform_control == T) {
        return(list("edge_center_dt"=combined_edge_center_dt,"control_dist_dt"=combined_control_dist_dt,"edge_malignant_dist_dt"=combined_edge_malignant_dist_dt,"merged_pca_loadings"=merged_pca_loadings_dt,"normal_pca_loadings"=normal_pca$loadings,"combined_pca_loadings"=combined_pca$loadings, "normal_pca"=normal_pca$pca,"combined_pca"=combined_pca$pca) )
    } else {
        return( list("edge_center_dt"=combined_edge_center_dt,"edge_malignant_dist_dt"=combined_edge_malignant_dist_dt,"merged_pca_loadings"=merged_pca_loadings_dt,"normal_pca_loadings"=normal_pca$loadings,"combined_pca_loadings"=combined_pca$loadings, "normal_pca"=normal_pca$pca,"combined_pca"=combined_pca$pca) )  }
}

perform_pca_shuffle <- function(normal_cells,normal_pca_df,num_control_shuffles,sample_wise_cell_name_list=NULL) {
    control_dist_matrix <- matrix(0,nrow=length(normal_cells),ncol=num_control_shuffles)
    for (i in 1:num_control_shuffles) {
        if (i%%10 == 0) {
            print(paste0(i,"/",num_control_shuffles))
            flush.console()
        }
        if (is.null(sample_wise_cell_name_list)) { 
            shuffled_pca_df <- t(apply(normal_pca_df,1,sample))
        } else {
            shuffled_pca_df <- normal_pca_df
            for (sample_ in names(sample_wise_cell_name_list)) { 
                cells <- sample_wise_cell_name_list[[sample_]]         

                if (length(cells) == 1) {
                    shuffled_pca_df[cells,] <- normal_pca_df[cells,]
                } else {
                    shuffled_pca_df[cells,] <- t(apply(normal_pca_df[cells,],1,sample))
                }
            }

        }
        ret <- self_pca_distances( normal_cells, shuffled_pca_df )
        control_dist_matrix[,i] = ret$dist_from_medoid
    }
    colnames(control_dist_matrix) <- paste("control",1:num_control_shuffles,sep="_")
    control_dist_dt <- data.table( control_dist_matrix )

    return(control_dist_dt)
}

perform_edge_center_shuffle <- function(edge_center_dt,normal_cell_type,malignant_cell_type,num_control_shuffles=100) {
    edge_center_dist_ratio <- mean(edge_center_dt[cell_category == "edge",dist_from_malignant_medoid])/mean(edge_center_dt[cell_category == "center",dist_from_malignant_medoid])
    edge_malignant_dist_dt <- data.table( normal_cell_type = normal_cell_type, malignant_cell_type = malignant_cell_type, edge_center_dist_ratio=edge_center_dist_ratio, dist_type="actual" )

    temp_dt <- copy(edge_center_dt[,.(cell_category, dist_from_malignant_medoid)])
    control_edge_center_dist_ratios <- vector(length=num_control_shuffles)
    for (shuffle_num in 1:num_control_shuffles) {
        temp_dt$cell_category <- sample(temp_dt$cell_category)
        control_edge_center_dist_ratios[shuffle_num] <- mean(temp_dt[cell_category == "edge",dist_from_malignant_medoid])/mean(temp_dt[cell_category == "center",dist_from_malignant_medoid])
    }
    control_malignant_dist_dt <- data.table( edge_center_dist_ratio=control_edge_center_dist_ratios, dist_type="control" )
    control_malignant_dist_dt[,`:=`(normal_cell_type=normal_cell_type,malignant_cell_type=malignant_cell_type)]
    edge_malignant_dist_dt <- rbind(edge_malignant_dist_dt,control_malignant_dist_dt)

    return(edge_malignant_dist_dt)
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
                sparse_cell <- as(count_mat[,cell_name],"sparseMatrix" )
                rownames(sparse_cell) <- all_genes
                colnames(sparse_cell) <- cell_name
                resampled_mat <- RowMergeSparseMatrices( resampled_mat, sparse_cell )
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

resample_edge_cells <- function( seurat_obj, pdac_anno_dt, edge_center_dt, normal_cell_types, malignant_cell_types, feature_to_subsample, no_tumour_adjacent=T, recompute_edge=T, perform_control=T, num_pcs=50 ) {
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

            if (recompute_edge) {
                resampled_obj_meta_data_dt <- merge( resampled_obj_meta_data_dt, pair_edge_center_dt, by="cell.name", suffixes=c("_after_resampling","_before_resampling"))

                resampled_edge_info <- add_edge_center_annotation( resampled_seurat_obj, pdac_anno_dt, malignant_cell_types=c(malignant_cell_type_),
                    normal_cell_types=c(normal_cell_type_), num_pcs=num_pcs, perform_control=perform_control, no_tumour_adjacent=no_tumour_adjacent )

                resampled_obj_meta_data_dt <- merge( resampled_obj_meta_data_dt, resampled_edge_info$edge_center_dt, by="cell.name",suffixes=c("_after_resampling","_before_resampling")  )
                old_edge_center_dt <- copy(pair_edge_center_dt)[,.(cell.name,cell_category,dist_from_malignant_medoid)]
                new_edge_center_dt <- merge( old_edge_center_dt, resampled_edge_info$edge_center_dt[,.(cell.name,dist_from_malignant_medoid)], by="cell.name" )
                new_edge_center_dt[,`:=`(dist_from_malignant_medoid=dist_from_malignant_medoid.y,dist_from_malignant_medoid.x=NULL)]
                edge_malignant_dist_original_edge_dt <- perform_edge_center_shuffle(new_edge_center_dt,normal_cell_type_,malignant_cell_type_)
                resampled_edge_info$edge_malignant_dist_dt <- edge_malignant_dist_original_edge_dt
            }
        }
    }

    if (recompute_edge) {
        return(list("resampled_meta_data_dt"=resampled_obj_meta_data_dt,"edge_info"=resampled_edge_info,"seurat_obj"=resampled_seurat_obj))
    } else {
        return(list("resampled_meta_data_dt"=resampled_obj_meta_data_dt,"seurat_obj"=resampled_seurat_obj))
    }
}

process_progenitor_data <- function() {
    progenitor_count_mat <- read_gene_exp_mat(file.path(base_path,"GSM4194789_TMM_counts_CPM.csv.gz"))

    progenitor_seurat_obj <- CreateSeuratObject( progenitor_count_mat )

    progenitor_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(progenitor_seurat_obj, pattern = "^MT-")
    dying_cells <- WhichCells( progenitor_seurat_obj, expression='percent.mt > 10' )
    
    #Seurat boilerplate routines to get UMAP coordinates from a read count matrix.
    progenitor_seurat_obj <- NormalizeData(progenitor_seurat_obj) 
    progenitor_seurat_obj <- ScaleData(progenitor_seurat_obj) %>% FindVariableFeatures(.,nfeatures = 1000) %>% RunPCA(.,dims=50) %>% RunUMAP(.,dims=1:50) 

    #Boilerplate code for running doubletFinder.
    homotypic.prop <- modelHomotypic(progenitor_seurat_obj@meta.data$cluster)
    nExp_poi <- round(0.05*nrow(progenitor_seurat_obj@meta.data))

    sample_param_sweep <- paramSweep_v3(progenitor_seurat_obj, PCs = 1:50, sct = FALSE)
    sample_param_sweep_summary <- summarizeSweep(sample_param_sweep, GT = FALSE)
    sweep_dt <- data.table( find.pK( sample_param_sweep_summary ) )
    optimal_pK <- as.double(as.vector(sweep_dt[order(-BCmetric),][1]$pK))
    progenitor_seurat_obj <- doubletFinder_v3(progenitor_seurat_obj, 
                               PCs = 1:50, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, 
                                          reuse.pANN = FALSE, sct = FALSE)

    sample_meta_data_dt <- data.table( progenitor_seurat_obj@meta.data, keep.rownames = T ) %>% setnames(.,"rn","cell.name")
    doublet_class_col <- paste("DF.classifications_0.25",optimal_pK,nExp_poi,sep="_")
    doublet_info_dt <- sample_meta_data_dt[,c("cell.name",doublet_class_col),with=F] %>%  setnames(.,doublet_class_col,"doublet_class")

    cell_doublets <- progenitor_meta_data_dt[cell.name %in% doublet_info_dt[doublet_class == "Doublet",cell.name],cell.name]
    
    cells_to_retain <- setdiff( Cells(progenitor_seurat_obj), c(dying_cells,doublet_cells) )
    #Removing doublets and dying cells, and re-running the standard Seurat pipeline. 
    progenitor_seurat_obj <- subset( progenitor_seurat_obj, cells=cells_to_retain )
    progenitor_meta_data_dt <- data.table( progenitor_seurat_obj@meta.data, keep.rownames=T ) %>% setnames(.,"rn","cell.name")

    progenitor_seurat_obj <- NormalizeData(progenitor_seurat_obj) %>% ScaleData(.) %>% 
    FindVariableFeatures(.,nfeatures = 2000) %>% RunPCA(.,dims=50) %>% FindNeighbors( ., reduction = "pca", k.param=30 ) %>%
    FindClusters(.,resolution=0.1) %>% RunUMAP(.,dims=1:50)
    progenitor_meta_data_dt <- data.table( progenitor_seurat_obj@meta.data, keep.rownames=T ) %>% setnames(.,"rn","cell.name")
    fwrite( progenitor_meta_data_dt, file.path(base_path,"progenitor_annotation.tsv"), sep="\t", row.names=F,quote=F)

    return(progenitor_seurat_obj)
}
 
compute_edge_frac <- function(age_data) {
    info_tbl <- age_data[,.N,by=list(age,cell_type)] %>% 
    tidyr::pivot_wider(names_from=c("cell_type"),names_prefix="N_",values_from=c("N")) %>% 
    tidyr::replace_na(list("N_edge"=0)) %>% mutate(frac_edge=N_edge/(N_edge+N_non_edge)) %>% dplyr::select(age,frac_edge) %>% arrange(age)
    
    return(info_tbl)
}

gene_assign <- function( regulatory_dt, assign_type="overlap", dist_threshold=8000 ) {
    granges_obj <- makeGRangesFromDataFrame( regulatory_dt, keep.extra.columns=T)

    promoters_txdb <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene,upstream=dist_threshold,downstream=400) %>% trim(.)

    mcol_df <- as.data.frame(mcols(promoters_txdb))
    mcol_df$gene_type <- "Expressed"
    mcols(promoters_txdb) <- mcol_df
    rm(mcol_df) 

    promoter_dt <- data.table( as.data.frame(promoters_txdb))
    
    overlap_obj <- findOverlaps( granges_obj, promoters_txdb )
    regulatory_dt[queryHits(overlap_obj),`:=`(TXNAME=promoter_dt[subjectHits(overlap_obj),tx_name])]

    regulatory_dt <- merge( regulatory_dt, cds_tx_dt[,.(TXNAME,gene_type,SYMBOL)], by="TXNAME", all.x=T )[,!c("TXNAME")]
    return(regulatory_dt)
}

compute_gene_set_AUCell_scores <- function(gene_exp_mat,gene_set,age_info_dt,auc_thr=NULL,threshold_type="Global_k1") {
    #gene_sets <- list("edge"=acinar_edge_genes)
    gene_set_name <- names(gene_set)
    cells_rankings <- AUCell_buildRankings(gene_exp_mat, nCores=6, plotStats=F)
    cells_AUC <- AUCell_calcAUC(gene_set, cells_rankings)
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

    aucell_scores_mat <- t(getAUC(cells_AUC))
    aucell_dt <- data.table(aucell_scores_mat,keep.rownames=T) %>% setnames("rn","cell.name") %>% #value=paste(gene_set_name,"AUC",sep="_") ) %>%
     merge( age_info_dt[,.(cell.name,age)], by="cell.name")
    if (is.null(auc_thr)) {
        auc_thr <- cells_assignment[[gene_set_name]]$aucThr$thresholds[threshold_type,"threshold"]
    }

    aucell_dt[,cell_type:=ifelse(get(gene_set_name) > auc_thr,gene_set_name,paste("non",gene_set_name,sep="_"))]

    #info_tbl <- compute_edge_frac( aucell_dt )#[,.N,by=c("age","cell_type")] )
    return(list("aucell_dt"=aucell_dt,"auc_threshold"=auc_thr))
}

load_pathways <- function() {
    pathways_h <- gmtPathways(file.path(base_path,"h.all.v7.0.symbols.gmt"))

    cancersea_gene_sets <- list()
    signature_file_paths <- list.files(file.path(base_path,"cancer-sea"),full.names=T)
    for (idx in 1:length(signature_file_paths)) {
        signature_genes_dt <- fread(signature_file_paths[idx])
        temp <- tail(unlist(strsplit(signature_file_paths[idx],"/")),1)
        signature_name <- paste( "CancerSEA", unlist( strsplit( temp, "\\." ) )[[1]], sep="_" )
        signature_gene_names <- signature_genes_dt$GeneName
        cancersea_gene_sets[[signature_name]] <- signature_gene_names
    }
    pathways_all <- c(pathways_h,cancersea_gene_sets)

    num_pathways <- length(pathways_all)

    for (pathway in names(pathways_all)) {
        if (length(pathways_all[[pathway]]) == 0) {
            pathways_all[[pathway]] <- NULL
        }
    }
    return(pathways_all)
}

compute_enrichment <- function(foreground_genes,all_genes,background_genes=NULL,pathways=NULL) { 
    if (is.null(pathways)) { 
        pathways <- load_pathways()
        pathways <- lapply( pathways, function(pathway_genes) {return(pathway_genes[pathway_genes %in% all_genes])})
    }
    fisher_enrichment_dt <- data.table(pathway=names(pathways),p_value=-1,odds_ratio=-1)
    if (is.null(background_genes)) {
        background_genes <- setdiff(all_genes,foreground_genes)
        #background_genes <- unique(unlist(pathways))
    }

    for (pathway in names(pathways)) {
        pathway_genes <- pathways[[pathway]]
        non_pathway_genes <- setdiff(all_genes,pathway_genes)

        num_in_pathway_and_foreground <- intersect(pathway_genes,foreground_genes) %>% length
        num_in_pathway_and_not_foreground <- intersect(pathway_genes,background_genes) %>% length
        num_not_in_pathway_and_foreground <- intersect(non_pathway_genes,foreground_genes) %>% length
        num_not_in_pathway_and_not_foreground <- intersect(non_pathway_genes,background_genes) %>% length
        fisher_mat <- matrix(c(num_in_pathway_and_foreground,num_in_pathway_and_not_foreground,
                              num_not_in_pathway_and_foreground,num_not_in_pathway_and_not_foreground),
                             nrow=2,ncol=2,byrow=T)
        
        test_res <- fisher.test(fisher_mat,alternative="g")
        pathway_ <- pathway
        fisher_enrichment_dt[pathway==pathway_,`:=`(p_value=test_res$p.value, odds_ratio=test_res$estimate,
        num_p_fg=num_in_pathway_and_foreground,num_p_bg=num_in_pathway_and_not_foreground,
        num_not_p_fg=num_not_in_pathway_and_foreground,num_not_p_bg=num_not_in_pathway_and_not_foreground)]
    }
    fisher_enrichment_dt[,q_value:=p.adjust(p_value)]

    return(fisher_enrichment_dt)
}

get_edge_signatures <- function(subset_seurat_obj, vars_to_regress=NULL,filter=NULL) {
    DefaultAssay(subset_seurat_obj) <- "RNA"
    gene_exp_mat <- subset_seurat_obj[["RNA"]]@data
    cell_meta_data_df <- subset_seurat_obj@meta.data
    gene_meta_data_df <- data.frame(gene_name=rownames(subset_seurat_obj))
    mast_obj <- MAST::FromMatrix(as.matrix(gene_exp_mat),cData=cell_meta_data_df,
                                fData=gene_meta_data_df)

    if (is.null(vars_to_regress)) { 
        mast_results_obj <- MAST::zlm(~cell_category, mast_obj)
    } else {
        f <- as.formula( paste("~ cell_category",paste(vars_to_regress,collapse="+"),sep="+" ))
        mast_results_obj <- MAST::zlm(f, mast_obj)
    }

    cond_obj <- summary(mast_results_obj,doLRT="cell_categoryedge")

    with_batch_p_value_dt <- cond_obj$datatable %>% dplyr::filter(component == "H" & contrast == "cell_categoryedge") %>% 
    mutate(fdr=p.adjust(`Pr(>Chisq)`))
    with_batch_fc_dt <- cond_obj$datatable %>% dplyr::filter(component == "logFC" & contrast == "cell_categoryedge") %>% arrange(-coef)

    return(list("fc_dt"=with_batch_fc_dt,"p_value_dt"=with_batch_p_value_dt))
}

load_liver_data <- function() {
    ident_vec <- c("1"="Hep1","3"="Hep2","6"="Hep4","14"="Hep5","15"="Hep6","5"="Hep3","17"="Cholangiocytes")
    normal_liver_seurat_obj <- readRDS(file.path(base_path,"HumanLiver_NameYourCellSubset.rds")) %>% NormalizeData %>%
    subset(subset=res.0.8 %in% c(1,3,6,14,15,5,17)) %>% RenameIdents(ident_vec)
    meta_data_df <- merge( tibble::rownames_to_column(normal_liver_seurat_obj@meta.data,"cell.name"),
                          tibble::enframe(ident_vec,name="res.0.8",value="Type"), by="res.0.8" ) %>%
    tibble::column_to_rownames("cell.name")
    normal_liver_seurat_obj <- AddMetaData(normal_liver_seurat_obj,meta_data_df)

    hcc_mat_1 <- Read10X(file.path(base_path,"GSE125449_Set1"))
    set_1_meta_data_dt <- fread(file.path(base_path,"GSE125449_Set1","GSE125449_Set1_samples.txt.gz"))
    hcc_mat_2 <- Read10X(file.path(base_path,"GSE125449_Set2"))
    set_2_meta_data_dt <- fread(file.path(base_path,"GSE125449_Set2","GSE125449_Set2_samples.txt.gz"))
    set_2_meta_data_dt <- fread(file.path(base_path,"GSE125449_Set2","GSE125449_Set2_samples.txt.gz"))

    hcc_obj_1 <- CreateSeuratObject( hcc_mat_1, 
                                    meta.data = tibble::column_to_rownames(set_1_meta_data_dt,"Cell Barcode"),
                                   project="Set1")
    hcc_obj_2 <- CreateSeuratObject( hcc_mat_2, meta.data = tibble::column_to_rownames(set_2_meta_data_dt,"Cell Barcode"),
                                   project="Set2")

    merged_liver_obj <- merge( x=hcc_obj_1, y=list(hcc_obj_2,normal_liver_seurat_obj ))
    rm(hcc_obj_1)
    rm(hcc_obj_2)
    rm(normal_liver_seurat_obj)
    merged_liver_anno_dt <- data.table( merged_liver_obj@meta.data, keep.rownames=T) %>% setnames("rn","cell.name")

    return(list("seurat_obj"=merged_liver_obj,"anno_dt"=merged_liver_anno_dt))
}

load_lung_data <- function() {
    luad_anno_dt <- fread("../data/GSE131907_Lung_Cancer_cell_annotation.txt.gz")
    sub_cell_types <- c("AT1","AT2","tS2","tS1","tS3","Ciliated","Club")
    #sub_cell_types <- c("AT1","AT2","MAST","NK","Treg","Monocytes","COL14A1+ matrix FBs","Ciliated","Club")
    all_samples <- luad_anno_dt$Sample %>% unique
    luad_samples_to_use <- all_samples[grepl("LUNG_N",all_samples) | grepl("LUNG_T",all_samples)]
    sub_luad_anno_dt <- luad_anno_dt[Sample %in% luad_samples_to_use & Cell_subtype %in% sub_cell_types,]
    cells_to_use <- sub_luad_anno_dt$Index

    luad_gene_exp_dt <- fread("../data/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
                              select=c("Index",cells_to_use))
    luad_gene_exp_mat <- luad_gene_exp_dt[,!c("Index")] %>% as.matrix
    rownames(luad_gene_exp_mat) <- luad_gene_exp_dt$Index
    rm(luad_gene_exp_dt)
    sub_luad_anno_dt[,cluster:=Cell_subtype]#fifelse(Cell_subtype %in% c("tS1","tS2","tS3"),"malignant",Cell_subtype)]
    luad_seurat_obj <- CreateSeuratObject( luad_gene_exp_mat, 
                                          meta.data=tibble::column_to_rownames(sub_luad_anno_dt,"Index"))
    rm(luad_gene_exp_mat)
    sub_luad_anno_dt[,sample_type:=fifelse(Sample_Origin == "nLung","normal","tumor")]
    setnames(sub_luad_anno_dt,"Index","cell.name")
    luad_seurat_obj <- SetIdent(luad_seurat_obj,value="cluster")

    return(list("seurat_obj"=luad_seurat_obj,"anno_dt"=sub_luad_anno_dt))
}

get_pca_correlations <- function(normal_pca,combined_pca) {
    pca_dt <- merge( data.table(combined_pca,keep.rownames=T), 
                    data.table(normal_pca,keep.rownames=T), by="rn") 

    pca_cor_mat <- matrix(0,nrow=ncol(normal_pca),ncol=ncol(combined_pca),
                         dimnames=list(colnames(normal_pca),colnames(combined_pca)))
    p_value_mat <- copy(pca_cor_mat)

    for (row_idx in 1:ncol(normal_pca)) {
        for (col_idx in 1:ncol(combined_pca)) {
            cor_info <- cor.test( pca_dt[[paste("normal_PC",row_idx,sep="_")]],pca_dt[[paste("combined_PC",col_idx,sep="_")]], method="pearson")
            pca_cor_mat[row_idx,col_idx] <- cor_info$estimate
            p_value_mat[row_idx,col_idx] <- cor_info$p.value
        }
    }

    pca_cor_df <- merge(reshape2::melt(pca_cor_mat,value.name="cor"),
                        reshape2::melt(p_value_mat,value.name="p_value"), by=c("Var1","Var2") ) %>% 
    mutate(q_value=p.adjust(p_value))
    return(pca_cor_df)
}

run_edge_pipeline <- function( seurat_obj, malignant_cell_types=NULL, normal_cell_types=NULL, num_pcs=50, perform_control = T, num_control_shuffles = 100, no_tumour_adjacent=F, num_var_features=1000, assay_to_use="RNA", ident_to_use="cluster", edge_cell_fraction=0.1, pipeline_variant="classic",features_prefix=NULL,feature_value_thresh=0.1,sample_info_column=NULL) {
    malignant_cell_types <- c(malignant_cell_types)
    normal_cell_types <- c(normal_cell_types)
   if (pipeline_variant == "classic") {
       edge_info <- add_edge_center_annotation_classic( seurat_obj, malignant_cell_types=malignant_cell_types,
       normal_cell_types=normal_cell_types, num_pcs=num_pcs, perform_control=perform_control,
       no_tumour_adjacent=no_tumour_adjacent, num_var_features=num_var_features,
       assay_to_use="RNA",ident_to_use=ident_to_use,edge_cell_fraction=edge_cell_fraction,
       sample_info_column=sample_info_column)
   } else if (pipeline_variant == "feature-selection") {
        print("Feature selection mode")
        features_file_name <- paste(normal_cell_types,malignant_cell_types,assay_to_use,"analysis.rds",sep="_")
        if (!is.null(features_prefix)) {
            features_file_name <- paste(features_prefix,features_file_name,sep="_")
        }

        if (!is.null(sample_info_column)) { 
            features_file_name <- paste("sample_aware",features_file_name,sep="_")
        }

        if (length(normal_cell_types) > 1 || length(malignant_cell_types) > 1) {
            edge_info_list <- list()
        }
        for (normal_cell_type in normal_cell_types) {
            for (malignant_cell_type in malignant_cell_types) {
                if (file.exists(features_file_name)) {
                   edge_features <- readRDS(features_file_name) 
                } else {
                    edge_features <- select_edge_center_features( seurat_obj,
                    malignant_cell_types=c(malignant_cell_type), normal_cell_types=c(normal_cell_type), num_pcs=num_pcs,
                    perform_control=perform_control, no_tumour_adjacent=no_tumour_adjacent,
                    num_var_features=num_var_features,
                    assay_to_use="RNA",ident_to_use=ident_to_use,edge_cell_fraction=edge_cell_fraction,
                    sample_info_column=sample_info_column)
                    saveRDS(edge_features,features_file_name)
                }

                skewness_pc <- edge_features$skewness_p_value[match(paste("normal_PC",1:num_pcs,sep="_"),Var1),][q_value < feature_value_thresh & normal_cell_type == eval(normal_cell_type),Var1] %>% head(1)
                if (length(skewness_pc) == 0)
                    return(NULL)

                edge_distance_pc <- edge_features$edge_center_p_value[q_value < feature_value_thresh & normal_cell_type == eval(normal_cell_type) & normal_pc == skewness_pc,pc] %>% head(1)

                if (length(edge_distance_pc) == 0)
                    return(NULL)

               edge_info <- add_edge_center_annotation_classic( seurat_obj, malignant_cell_types=c(malignant_cell_type),
               normal_cell_types=c(normal_cell_type), num_pcs=num_pcs, perform_control=perform_control,
               normal_pcs_to_use=skewness_pc, pooled_pcs_to_use=edge_distance_pc,
               no_tumour_adjacent=no_tumour_adjacent, num_var_features=num_var_features,
               assay_to_use="RNA",ident_to_use=ident_to_use,edge_cell_fraction=edge_cell_fraction)

                pca_cor_df <- get_pca_correlations( edge_info$normal_pca, edge_info$combined_pca ) %>% mutate(q_value=p.adjust(p_value)) %>% arrange(q_value)
                edge_info[["pca_cor_df"]] <- pca_cor_df
                if (length(normal_cell_types) > 1 || length(malignant_cell_types) > 1) {
                    edge_info_list[[paste(normal_cell_type,malignant_cell_type,sep="-")]] <- edge_info
                }   
            }
        }
    }

    if (length(normal_cell_types) > 1 || length(malignant_cell_types) > 1) {
        return(edge_info_list) 
    } else {
        return(edge_info)
    }
}

gene_exp_regression <- function(gene_exp_mat,to_regress_df,regressors,coefs_to_return) {
    glm_formula <- as.formula( paste("z_score", paste(regressors,collapse="+"), sep="~"))
    coefs_to_return <- c(coefs_to_return)
    p_value_vars <- paste(coefs_to_return,"p_value",sep="_")
    coef_vars <- paste(coefs_to_return,"coef",sep="_")
    age_cor_dt <- data.table(gene = rownames(gene_exp_mat) )
    age_cor_dt[,eval(p_value_vars):=-1]
    age_cor_dt[,eval(coef_vars):=-1]

    to_regress_df <- tibble::rownames_to_column(to_regress_df,"name")
    for (gene_name in rownames(gene_exp_mat)) {
        gene_exp_vec <- gene_exp_mat[gene_name,]
        if (sum(gene_exp_vec > 0) == 0 )
            next
        gene_exp_df <- tibble::enframe(gene_exp_vec,value="z_score") %>% merge( ., to_regress_df, by="name" )
        glm_obj <- glm(glm_formula,data=gene_exp_df)
        coef_table <- coef(summary(glm_obj))
        for (idx in 1:length(coefs_to_return)) {
            glm_p_value <- coef(summary(glm_obj))[,4][coefs_to_return[idx]]
            estimate <- coef(summary(glm_obj))[,1][coefs_to_return[idx]]
            age_cor_dt[gene == gene_name,eval(coef_vars[idx]):=estimate]
            age_cor_dt[gene == gene_name,eval(p_value_vars[idx]):=glm_p_value]
        }
    }
    q_value_vars <- paste(coefs_to_return,"q_value",sep="_")
    for (idx in 1:length(coefs_to_return)) {
        age_cor_dt <- age_cor_dt[get(p_value_vars[idx]) != -1,]
        q_values <- p.adjust(age_cor_dt[[p_value_vars[idx]]])
        age_cor_dt[,eval(q_value_vars[idx]):=q_values]
    }

    return(age_cor_dt)
}

load_GSE81547_mutations <- function(aucell_dt_) { 
    mutations_dt <- fread(file.path(base_path,"binary.SC"))
    sra_info_df <- fread(file.path(base_path,"sra_info_table.csv"))
    merged_edge_sra_dt <- merge( aucell_dt_, sra_info_df %>% dplyr::select(Run,`Sample Name`), by.x="cell.name", 
                                by.y="Sample Name") %>% mutate(age=as.character(age)) %>% data.table
    merged_edge_sra_dt$age <- factor( merged_edge_sra_dt$age, as.character(sort(as.integer(unique(merged_edge_sra_dt$age)),decreasing=F)))

    mutation_names <- setdiff( names(mutations_dt),c("cellIDxmutID", "age" ))
    mutation_coord_dt <- t(sapply( str_split(mutation_names,"\\."), unlist ))[,c(1,2,3,4,5)] %>% data.table
    mutation_coord_dt <- mutation_coord_dt[,.(gene=V1,chr=V2,start=as.integer(V3),ref=V4,alt=V5,dbsnp=F)]
    chromosomes <- c(paste0("chr",1:22),"X","Y","M")#unique(mutation_coord_dt$chr)
    for (chrom in chromosomes) {
        #file_path <- paste( paste0("../data/dbsnp_chrom_wise/",chrom), "vcf.gz",sep=".")
        file_path <- paste( file.path(base_path,"dbsnp_chrom_wise",chrom), "vcf.gz",sep=".")
        if (file.exists(file_path)) {
            vcf_dt <- fread(file_path,
                            header=F,select=c("V1","V2"),showProgress=T)
            start_coords <- mutation_coord_dt[chr == chrom,start]
            mutation_coord_dt[chr == chrom & start %in% vcf_dt$V2,dbsnp:=T]
            print(chrom)
            flush.console()
        }
    }
    mutations_to_keep <- mutation_coord_dt[dbsnp == F,] %>% dplyr::select(-dbsnp) %>% 
    mutate(mut_name=paste(gene,chr,start,ref,alt,sep=".")) %>% pull(mut_name)
    retained_mutations_dt <- mutations_dt[,c("cellIDxmutID",mutations_to_keep),with=F]
    mutations_per_cell_df <- retained_mutations_dt %>% tidyr::pivot_longer(cols=all_of(mutations_to_keep),values_to="mut_status") %>% 
    dplyr::filter(mut_status == 1) %>% merge(merged_edge_sra_dt[,.(Run,cell_type,age)],by.x="cellIDxmutID",by.y="Run") %>% 
    dplyr::select(-c(mut_status))
    shared_mutations_across_donors <- mutations_per_cell_df %>% dplyr::select(-cellIDxmutID,-cell_type) %>% unique %>% arrange(name) %>%
    dplyr::group_by(name) %>% dplyr::count() %>% dplyr::filter(n>1) %>% pull(name)
    retained_mutations_dt <- retained_mutations_dt[,!shared_mutations_across_donors,with=F]
    mutations_per_cell_df <- mutations_per_cell_df %>% dplyr::filter(!name %in% shared_mutations_across_donors)
}

process_GSE81547 <- function() { 
    gene_exp_dir <- file.path(base_path,"GSE81547")
    out_path <- file.path(gene_exp_dir,"GSE81547_mat.tsv.gz")
    print(out_path)
    
    if (!file.exists(out_path)) {
        file_paths = list.files(gene_exp_dir,full.names=T)
        gene_exp_dt <- data.table()
        col_idx <- 1
        for (file_path in file_paths) {
            if (grepl(".gz",file_path)) {
                flush.console()
                dt <- fread( file_path )
                gsm_num <- str_match( file_path, "GSM[0-9][0-9]*" )[1]
                if (nrow(gene_exp_dt) == 0 ) {
                    gene_exp_dt <- dt %>% setnames(.,"V2",gsm_num)
                } else {
                    gene_exp_dt <- cbind( gene_exp_dt, dt[,list(V2)] ) %>% setnames(.,"V2",gsm_num)
                }
                col_idx <- col_idx + 1
            }
        }
        fwrite(gene_exp_dt,out_path,quote=F,row.names=F,sep="\t")
    }

    gene_exp_mat <- read_gene_exp_mat(out_path)
    aging_obj <- create_full_seurat_object( gene_exp_mat )
    aging_obj <- NormalizeData( aging_obj )

    rm(gene_exp_mat)

    geo_info <- getGEO("GSE81547",getGPL=F)
    geo_info_df <- pData(geo_info$GSE81547_series_matrix.txt.gz)
    donor_age <- str_match(geo_info_df[,"title"],"^[0-9][0-9]*")
    info_to_add_df <- data.frame( age=as.integer(donor_age), row.names=Cells(aging_obj))

    anno_dt <- fread(file.path(base_path,"acinar_scrna_seq.txt"))
    cell_types <- anno_dt[match(`Sample Name`,rownames(info_to_add_df)),inferred_cell_type]
    info_to_add_df$cell_type <- cell_types

    aging_obj <- AddMetaData( aging_obj, info_to_add_df )

    return(aging_obj)
}

process_GSE85241 <- function() {
    gene_exp_mat <- read_gene_exp_mat(file.path(base_path,"GSE85241_cellsystems_dataset_4donors_updated.csv.gz"))

    processed_gene_names <- gsub("__.*","",rownames(gene_exp_mat))
    names(processed_gene_names) <- rownames(gene_exp_mat)
    duplicate_gene_names <- processed_gene_names[duplicated( processed_gene_names ) == T]
    unique_gene_names <- processed_gene_names[!processed_gene_names %in% duplicate_gene_names]
    new_gene_exp_mat <- gene_exp_mat[names(unique_gene_names),]
    rownames(new_gene_exp_mat) <- unique_gene_names
    new_mat <- matrix(0,nrow=length(duplicate_gene_names),ncol=ncol(new_gene_exp_mat),
                     dimnames=list(duplicate_gene_names,colnames(new_gene_exp_mat)))
    for (gene in duplicate_gene_names) {
        entries <- rownames(gene_exp_mat)[grepl(gene,rownames(gene_exp_mat))]
        new_mat[gene,] <- rowSums(gene_exp_mat[entries,])
    }
    new_gene_exp_mat <- rbind(new_gene_exp_mat,new_mat)
    cel_seq_obj <- create_full_seurat_object( new_gene_exp_mat )

   cel_seq_obj <- NormalizeData(cel_seq_obj) %>% ScaleData(.) %>% 
    FindVariableFeatures(.,nfeatures = 2000) %>% RunPCA(.,dims=50) %>% FindNeighbors( ., reduction = "pca", k.param=30 ) %>%
    FindClusters(.) %>% RunUMAP(.,dims=1:50)
    cel_seq_meta_data_dt <- data.table( cel_seq_obj@meta.data, keep.rownames=T ) %>% setnames(.,"rn","cell.name")
    prss1_dt <- FetchData( cel_seq_obj, vars=c("seurat_clusters","PRSS1")) %>% data.table(.)
    prss1_dt[,mean_prss1:=mean(PRSS1),by="seurat_clusters"]
    prss1_dt <- prss1_dt[order(-mean_prss1),.(seurat_clusters,mean_prss1)] %>% unique(.)
    acinar_cluster <- prss1_dt[1,seurat_clusters]
    acinar_cells <- cel_seq_meta_data_dt[seurat_clusters == prss1_dt[1,seurat_clusters],cell.name]

    age_meta_data_dt <- data.table( donor_id=c("D28","D29","D30","D31","D16","D25"),
                                  age=c(54,23,48,59,53,30),
                                  sex=c("M","M","F","M","M","M"))

    cel_seq_meta_data_dt$donor_id <- gsub("-[0-9][0-9]*","",cel_seq_meta_data_dt$orig.ident)

    cel_seq_meta_data_dt <- merge( cel_seq_meta_data_dt, age_meta_data_dt, by="donor_id")

    return(list("seurat_obj"=cel_seq_obj,"meta_data"=cel_seq_meta_data_dt,"acinar_cells"=acinar_cells))
 
}

process_GSE85241_mutations <- function(vano_aucell_dt_) {
    mutations_dt <- fread(file.path(base_path,"GSE85241.csv")) %>% setnames("V1","cell.name")
    sra_info_dt <- fread(file.path(base_path,"GSE85241_Sra.txt")) %>% dplyr::select(Run,`Sample Name`)
    gse85241 <- getGEO("GSE85241")
    geo_info_df <- pData(gse85241$GSE85241_series_matrix.txt.gz) %>% dplyr::select(orig.ident=description,`Sample Name`=geo_accession)

    mutations_dt <- mutations_dt %>%
     mutate(Run=as.vector(unlist(str_match(cell.name,"SRR[0-9][0-9]*"))),
           cell_idx=as.vector(unlist(str_match(cell.name,"[0-9][0-9]*$"))) %>%
           gsub("^0*","",.)) %>% #%>% dplyr::select(cell.name,Run,cell_idx) %>%
    merge(.,sra_info_dt,by="Run") %>% merge(.,geo_info_df,by="Sample Name") %>% 
    mutate(cell.name=paste(orig.ident,cell_idx,sep="_")) %>%
    dplyr::select(-all_of(c("Run","Sample Name","cell_idx","orig.ident")))

    mutations <- colnames(mutations_dt) %>% setdiff(.,"cell.name")
    mutations_per_cell_df <- merge(cel_seq_meta_data_dt[cell.name %in% vano_acinar_cells,.(age,cell.name)],
    vano_aucell_dt_[,.(cell.name,cell_type)], by="cell.name") %>%
    merge(.,mutations_dt,by="cell.name") %>% 
    tidyr::pivot_longer(.,all_of(mutations),names_to="mutation",values_to="mut_status") %>%
    dplyr::filter(mut_status %in% c(1,3)) %>% dplyr::select(-mut_status)

    mutations_to_keep <- mutations_per_cell_df %>% dplyr::select(age,mutation) %>% unique %>%
dplyr::group_by(mutation) %>% dplyr::count() %>% dplyr::filter(n==1) %>% pull(mutation)
    mutations_per_cell_df <- mutations_per_cell_df %>% dplyr::filter(mutation %in% mutations_to_keep)
    return(mutations_per_cell_df)
}

compute_pseudotime <- function( cds_obj, genes_to_use, shuffle_pca=F, random.seed=1 ) {
    cds_obj <- preprocess_cds(cds_obj, num_dim = 50, use_genes=genes_to_use,
                            scaling=T, norm_method="none", method="PCA",
                             alignment_group = "orig.ident",
                          residual_model_formula_str = "~S.Score + G2M.Score")

    if (shuffle_pca) {
        set.seed(random.seed)
        pca_mat <- reducedDims(cds_obj)$PCA
        cell_names <- rownames(pca_mat)
        shuffled_pca_mat <- apply(reducedDims(cds_obj)$PCA,2,sample)
        rownames(shuffled_pca_mat) <- cell_names
        reducedDims(cds_obj)$PCA <- shuffled_pca_mat
    }

    cds_obj <- reduce_dimension(cds_obj,reduction_method = "UMAP")
    cds_obj <- cluster_cells(cds_obj)
    cds_obj <- learn_graph( cds_obj, use_partition=F )

    return(cds_obj)
}

set_root <- function(cds_obj,root_candidates) {
   closest_vertex <- cds_obj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
   closest_vertex <- as.matrix(closest_vertex[colnames(cds_obj), ])
   root_pr_nodes <-  igraph::V(principal_graph(cds_obj)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[root_candidates,]))))]

   cds_obj <- order_cells( cds_obj, root_pr_nodes=root_pr_nodes)
    return(cds_obj)
}

compute_pseudotime_distance_ratio <- function( normal_df_, malignant_df_ ) {
    malignant_median <- min(malignant_df_$pseudotime)
    mean_edge_distance <- mean((normal_df_ %>% dplyr::filter(cell_category=="edge") %>% pull(pseudotime) - malignant_median)^2)
    mean_center_distance <- mean((normal_df_ %>% dplyr::filter(cell_category=="center") %>% pull(pseudotime) - malignant_median)^2)

    return(mean_edge_distance/mean_center_distance)
}

search_cosmic_database <- function(cosmic_all_dt,cosmic_census_dt,mutations) {
    mutation_status_dt <- data.table(mutation=mutations,mut_info="")
    prev_chrom <- ""

    #The short loop below sorts the mutations by chromosome in order to speed up the look-ups of mutations
    #in the large COSMIC table.
    chroms <- paste0("chr",c(1:22,"X","Y","M"))
    mutations_to_scan <- c()
    for (chrom in chroms) {
        muts <- mutations[grepl(paste0(chrom,"\\."),mutations)]
        mutations_to_scan <- c(mutations_to_scan,muts)
    }
    for (mutation_ in mutations_to_scan) {
        flush.console()
        mut_info <- str_split(mutation_,"\\.") %>% unlist
        chrom <- mut_info[2]
        pos <- mut_info[3]
        ref <- mut_info[4]
        alt <- mut_info[5]
        if (chrom != prev_chrom) {
            print(paste("Caching",chrom))
            chrom_dt <- cosmic_all_dt[`#CHROM` == chrom,]
            prev_chrom <- chrom
            flush.console()
        }
        cosmic_info <- chrom_dt[`#CHROM` == chrom & POS == pos & REF == ref &
                                    ALT == alt,INFO]
        if (length(cosmic_info) > 0) {
            mutation_status_dt[mutation == mutation_,mut_info:=c(cosmic_info)[1]]
        }
    }
    rm(chrom_dt)
    mutations_in_cosmic_dt <- mutation_status_dt[mut_info != "",]
    num_samples_per_mutation <- str_match(mutation_status_dt[mut_info != "",]$mut_info,"CNT=[0-9][0-9]*") %>%
    str_match(.,"[0-9][0-9]*") %>% as.vector
    mutations_in_cosmic_dt$num_samples <- as.numeric(num_samples_per_mutation)
    mutation_ids <- str_match(mutations_in_cosmic_dt$mut_info,"COSM[0-9][0-9]*") %>% as.vector

    print("Tier of mutations present in sample")
    census_results_dt <- cosmic_census_dt[LEGACY_MUTATION_ID %in% mutation_ids,]
    print(unique(census_results_dt[,MUTATION_SIGNIFICANCE_TIER]))

    return(list(mutations_in_cosmic_dt,census_results_dt))
}

run_pc_wise_edge_analysis <- function(seurat_obj,normal_cell_types,malignant_cell_type,ident_to_use,sample_info_column,
q_value_thresh=0.1) {
    edge_info_list <- list()
    pathways <- load_pathways()
    for (normal_cell_type in normal_cell_types) {
        print(normal_cell_type)
        edge_features <- readRDS( paste("sample_aware",normal_cell_type,malignant_cell_type,"RNA_analysis.rds",sep="_") )
        skewness_pcs <- edge_features$skewness_p_value[q_value < 0.1,Var1]
        edge_distance_pcs <- edge_features$edge_center_p_value[q_value < 0.1 & normal_pc %in% skewness_pcs,pc] %>% unique

        if (length(skewness_pcs) == 0 || length(edge_distance_pcs) == 0) {
            next
        }

        for (skewness_pc in skewness_pcs) {
            for (edge_distance_pc in edge_distance_pcs) {
                flush.console()
                edge_info <- add_edge_center_annotation_classic(seurat_obj,normal_cell_types=normal_cell_type,
                       malignant_cell_types=malignant_cell_type,ident_to_use=ident_to_use,
                                                                            normal_pcs_to_use=c(skewness_pc),
                                                                            pooled_pcs_to_use=c(edge_distance_pc),
                                        sample_info_column=sample_info_column, num_pcs=5)
                pca_cor_df <- get_pca_correlations( edge_info$normal_pca, edge_info$combined_pca ) %>%
                dplyr::filter(Var1 %in% skewness_pcs & Var2 %in% edge_distance_pcs) %>%  mutate(q_value=p.adjust(p_value))
                collinearity_q_value <- pca_cor_df %>% dplyr::filter(Var1 == skewness_pc & Var2 == edge_distance_pc) %>% pull(q_value)
                if (collinearity_q_value < 0.1) {
                    print(skewness_pc)
                    print(edge_distance_pc)
                    print("--------")
                    if (is.null(edge_info_list[[normal_cell_type]])) {
                        edge_info_list[[normal_cell_type]] <- list()
                    }

                    edge_info_list[[normal_cell_type]][[skewness_pc]] <- edge_info


                    logFC_file <- paste(normal_cell_type,paste(skewness_pc,"rds",sep="."),sep="_")
                    if (!file.exists(logFC_file)) {
                        subset_seurat_obj <- subset( seurat_obj, subset = !!sym(ident_to_use) == normal_cell_type ) %>%
    AddMetaData( edge_info$edge_center_dt[,.(cell.name,cell_category)] %>%
               tibble::column_to_rownames("cell.name") ) %>% SetIdent(value="cell_category") %>% NormalizeData %>%
                   CellCycleScoring(.,g2m.features = g2m.genes,s.features=s.genes )

                        if (length(unique(subset_seurat_obj@meta.data[[sample_info_column]])) == 1) {
                            ret <- get_edge_signatures( subset_seurat_obj, vars_to_regress=c("S.Score","G2M.Score"))
                        } else {
                            ret <-  get_edge_signatures( subset_seurat_obj, vars_to_regress=c("S.Score","G2M.Score",sample_info_column))
                        }
                        saveRDS( ret, logFC_file )
                    } else {
                        ret <- readRDS(logFC_file)
                    }

                    sorted_gene_vec <- ret$fc_dt[order(-coef),.(primerid,coef)] %>% tibble::deframe(.)
                    fgsea_dt <- fgsea( pathways, sorted_gene_vec, nperm=1000 )
                    edge_info_list[[normal_cell_type]][[skewness_pc]]$fgsea_dt <- fgsea_dt
                    edge_info_list[[normal_cell_type]][[skewness_pc]]$pca_cor_df <- pca_cor_df
                    break
                }

            }
        }
    }

    return(edge_info_list)
}

check_mutation_based_edge_cells <- function(mutations_per_cell_df,seurat_obj,aucell_dt,logfc_thresh=1,p_val_adj_thresh=0.1) {
    mutations <- unique(mutations_per_cell_df$mutation)
    all_cells <- unique(mutations_per_cell_df$cell.name)
    acinar_obj <- subset( seurat_obj, cells=all_cells )
    num_genes_above_thresh <- c()
    for (mutation_ in mutations) {
        with_mut_cells <- mutations_per_cell_df %>% dplyr::filter(mutation == mutation_) %>% pull(cell.name)
        without_mut_cells <- setdiff(all_cells,with_mut_cells)
        list_name <- paste(sort(with_mut_cells),collapse="_")
        
        df <- data.frame( cell.name=all_cells ) %>% 
        mutate(cell_type=ifelse(cell.name %in% with_mut_cells,"edge","non-edge")) %>% tibble::column_to_rownames("cell.name")
        
        acinar_obj <- AddMetaData( acinar_obj, df ) %>% SetIdent(value="cell_type")
        num_markers <- tryCatch( { FindMarkers( acinar_obj, ident.1="edge", logfc.thresh=logfc_thresh, min.cells.group = 1,
        test.use="LR",latent.vars=c("S.Score","G2M.Score","nFeature_RNA")) %>% 
        dplyr::filter(p_val_adj < p_val_adj_thresh) %>% nrow}, error=function(cond){return(0)} )
        num_genes_above_thresh <- c(num_genes_above_thresh,num_markers)
    }
    names(num_genes_above_thresh) <- mutations

    acinar_obj <- AddMetaData( acinar_obj, 
                              tibble::column_to_rownames(aucell_dt[,.(cell.name,cell_type)],"cell.name")) %>%
    SetIdent(value="cell_type")
    num_edge_markers <- FindMarkers( acinar_obj, ident.1="edge", logfc.thresh=logfc_thresh ) %>% 
    dplyr::filter(p_val_adj < 0.1) %>% nrow
    p_value <- 1 - pnorm( num_edge_markers, mean=mean(num_genes_above_thresh), sd=sd(num_genes_above_thresh))
    print("p-value of observed number of edge markers")
    print(p_value)


    return(list("num_edge"=num_edge_markers,"null_dist"=num_genes_above_thresh))
}

subsample_edge_non_edge_cells <- function(aucell_dt) {
    ages <- unique(aucell_dt$age)
    cells_to_retain <- c()
    for (age_ in ages) {
        edge_cells <- aucell_dt[age == age_ & cell_type == "edge",cell.name]
        non_edge_cells <- aucell_dt[age == age_ & cell_type == "non_edge",cell.name]
        num_cells_per_type <- aucell_dt[age == age_,.N,by=cell_type]
        num_edge <- num_cells_per_type[cell_type == "edge",N]
        num_non_edge <- num_cells_per_type[cell_type != "edge",N]
        if (length(num_non_edge) == 0) {
            num_non_edge = 0
        }
        
        if (length(num_edge) == 0) {
            num_edge = 0
        }

        if (num_non_edge > 0 && num_edge > 0){
            if (num_non_edge > num_edge) {
                non_edge_cells <- sample(aucell_dt[age == age_ & cell_type == "non_edge",cell.name],num_edge)
            } else if (num_edge > num_non_edge) {
                edge_cells <- sample(aucell_dt[age == age_ & cell_type == "edge",cell.name],num_non_edge)
            }
        }

        cells_to_retain <- c(cells_to_retain,c(edge_cells,non_edge_cells))
    }

    return(cells_to_retain)
}

