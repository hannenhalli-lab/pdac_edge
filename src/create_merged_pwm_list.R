jaspar2020_url <- "http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.zip"
jaspar_pfm_path <- file.path( "..", "data", "jaspar" )

download_jaspar_pfms <- function( jaspar_url ) {
    if (!dir.exists(jaspar_pfm_path)) {
        dir.create(jaspar_pfm_path)
    }
    jaspar_downloaded_file <- file.path( jaspar_pfm_path, "jaspar_pfms.zip")
    download.file( jaspar_url, destfile=jaspar_downloaded_file )
    unzip( jaspar_downloaded_file, exdir=jaspar_pfm_path )
    file.remove( jaspar_downloaded_file )
}

if (!dir.exists(jaspar_pfm_path))
    download_jaspar_pfms( jaspar2020_url )

hocomoco_pfm_path <- file.path( "..", "data", "hocomoco" )

download_hocomoco_pfms <- function( hocomoco_url ) {
    if (!dir.exists(hocomoco_pfm_path)) {
        dir.create(hocomoco_pfm_path)
    }
    hocomoco_downloaded_file <- file.path( hocomoco_pfm_path, "hocomoco_pfms.tar.gz")
    if (!file.exists(hocomoco_downloaded_file))
        download.file( hocomoco_url, destfile=hocomoco_downloaded_file )

    files_present <- list.files(hocomoco_pfm_path,full.names=T)
    if (length(files_present) == 0) {
        downloaded_files <- untar( hocomoco_downloaded_file, exdir=hocomoco_pfm_path, list=TRUE )
    } else {
        downloaded_files <- files_present
    }
    return(downloaded_files)
}
hocomocov11_url <- "https://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwm_HUMAN_mono.tar.gz"
hocomoco_pwm_paths <- download_hocomoco_pfms( hocomocov11_url )
hocomoco_pwm_paths <- file.path( hocomoco_pfm_path, hocomoco_pwm_paths )

hocomocov11_anno_url <- "https://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"
hocomoco_anno_path <- file.path( hocomoco_pfm_path, "annotation.tsv")
if (!file.exists(hocomoco_anno_path))
    download.file( hocomocov11_anno_url, hocomoco_anno_path )

hocomocov11_pvalue_thresh <- "https://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt"
hocomoco_p_value_path <- file.path( hocomoco_pfm_path, "pvalues.tsv")
if (!file.exists(hocomoco_p_value_path)) {
    download.file( hocomocov11_anno_url, hocomoco_p_value_path )
}
bg_nuc_freq<-c("A"=2.955e-01,"C"=2.045e-01,"G"=2.045e-01,"T"=2.955e-01)
create_hocomoco_pwm_list <- function( ) {
    anno_dt <- fread( hocomoco_anno_path )
    pwm_list <- PWMatrixList()

    for (pwm_path in hocomoco_pwm_paths) {
        if (!grepl(".pwm$",pwm_path))
            next

        motif_mat <- fread(pwm_path,skip=1) %>% as.matrix(.) %>% t(.)
        rownames(motif_mat) <- c("A","C","G","T")
        motif_name <- unlist( str_split( pwm_path, "/" ) ) %>% tail(.,n=1) %>% gsub(".pwm$","",.)
        tf_name <- anno_dt[Model == motif_name,`Transcription factor`]

        pwm_obj <- PWMatrix(ID=motif_name,name=tf_name,bg=bg_nuc_freq,profileMatrix=motif_mat)
        pwm_list[[tf_name]] <- pwm_obj
    }

    return(pwm_list)
}

merged_pwm_list <- create_hocomoco_pwm_list()

##Check if download was successful.
convert_jaspar_pfm_to_pwm <- function( ) {
    jaspar_pfm_paths <- list.files( jaspar_pfm_path, full.names=T )
    
    get_tf_name <- function( file_path ) {
        file_obj <- file(file_path,"r")
        jaspar_pfm_header <- readLines(file_obj,n=1)
        close(file_obj)
        header_chunks <- str_split_fixed( jaspar_pfm_header, "\t", n=2 )
        tf <- header_chunks[2] 
        motif_name <- header_chunks[1]
    
        return(list("tf"=tf,"motif_name"=motif_name))
    }
    
    for (jaspar_pfm_path in jaspar_pfm_paths) {
        ret <- get_tf_name( jaspar_pfm_path )
        tfs_present <- names(merged_pwm_list)

        if (ret$tf %in% tfs_present) {
            print(paste(ret$tf,"already present.Skipping."))
            flush.console()
            next
        }
        
        dt <- fread(jaspar_pfm_path,skip=1)
        last_col <- tail(names(dt),n=1)
        pfm_mat <- dt[,!c("V1","V2",last_col),with=F] %>% as.matrix(.)
        rownames(pfm_mat) <- c("A","C","G","T")
        pfm_obj <- PFMatrix(ID=ret$motif_name,name=ret$tf,bg=bg_nuc_freq,profileMatrix=pfm_mat)
        merged_pwm_list[[ret$tf]] <- toPWM(pfm_obj)
    }
    
    return(merged_pwm_list)
}

merged_pwm_list <- convert_jaspar_pfm_to_pwm()
