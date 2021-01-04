library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(universalmotif)
library(readr)

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
bg_nuc_freq<-c("A"=2.955e-01,"C"=2.045e-01,"G"=2.045e-01,"T"=2.955e-01)

##Check if download was successful.
convert_jaspar_pfm_to_pwm <- function( ) {
    jaspar_pfm_paths <- list.files( jaspar_pfm_path, full.names=T )
    merged_pwm_list <- list()
    
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

        dt <- read_table2(jaspar_pfm_path,skip=1,col_names=F,progress=F, col_types=cols()) %>% data.table(.)
        last_col <- tail(names(dt),n=1)
        pfm_mat <- dt[,!c("X1","X2",last_col),with=F] %>% as.matrix(.)
        rownames(pfm_mat) <- c("A","C","G","T")
        pfm_obj <- PFMatrix(ID=ret$motif_name,name=ret$tf,bg=bg_nuc_freq,profileMatrix=pfm_mat)
        merged_pwm_list[[ret$tf]] <- toPWM(pfm_obj)
    }
    
    return(merged_pwm_list)
}

merged_pwm_list <- convert_jaspar_pfm_to_pwm()

genome <- BSgenome.Hsapiens.UCSC.hg19
get_motif_score_thresh <- function( pwm_file_name, motif_p_value_thresh ) {
    ape_command <- paste("java -cp ape.jar ru.autosome.ape.FindThreshold", shQuote(pwm_file_name),
                         motif_p_value_thresh, "-b", paste( bg_nuc_freq, sep="," ), " | grep -v '^#'" )
    motif_score_thresh <- fread( cmd=ape_command )$V3

    return(motif_score_thresh)
}

create_motif_score_thresh_table <- function( pwm_list, motif_p_value_thresh ) {
    tf_names <- names(pwm_list)
    motif_p_value_dt <- data.table( tf=tf_names, score_thresh=-1 )
    score_thresh_vec <- c()
    for (tf_name in tf_names) {
        pwm_obj <- pwm_list[[tf_name]]
        print(tf_name)
        flush.console()
        pwm_dt <- data.table( t(Matrix(merged_pwm_list[[tf_name]])) )
        pwm_file_name <- paste(tf_name,"pwm",sep=".")
        fwrite( data.table( pwm_dt ), pwm_file_name, col.names=F, row.names=F, quote=F, sep="\t" )

        motif_score_thresh <- get_motif_score_thresh(pwm_file_name,motif_p_value_thresh)
        score_thresh_vec <- c(score_thresh_vec,motif_score_thresh)
        file.remove(pwm_file_name)
    }
    motif_p_value_dt[,score_thresh:=score_thresh_vec]

    return(motif_p_value_dt)
}

get_sarus_motif_scores <- function( file_name, pwm_file_name, motif_p_value_thresh ) {
#    if (!enhancer_search_mode) {
#        sarus_command <- paste("java -cp sarus-latest.jar ru.autosome.SARUS",
#                           shQuote(file_name), shQuote(pwm_file_name), "besthit", "--output-bed" )
#
#
#        dt <- fread( cmd=sarus_command )  #The 5th column, V5, is the motif score we need.
#        motif_scores <- dt$V5
#
#        return(motif_scores)
#    } else {
    motif_score_thresh <- get_motif_score_thresh( pwm_file_name, motif_p_value_thresh )
    sarus_command <- paste("java -cp sarus-latest.jar ru.autosome.SARUS",
                       shQuote(file_name), shQuote(pwm_file_name),  motif_score_thresh, "--output-bed" )

    dt <- fread( cmd=sarus_command )  #The 5th column, V5, is the motif score we need.
    if (nrow(dt) > 0 ) {
        sequences_containing_match <- unique(dt$V1)
        num_motif_matches_per_sequence <- rep(1,length(sequences_containing_match))
        names(num_motif_matches_per_sequence) <- sequences_containing_match
     } else {
        num_motif_matches_per_sequence <- NULL
    }
    return(num_motif_matches_per_sequence)
#    }
}

pwm_scan <- function( pwm_list, foreground_granges_obj, background_granges_obj=NULL,
                     do_enrichment=T, enhancer_search_mode=F, motif_p_value_thresh=NULL ) {
    set.seed(364002204)

    foreground_sequences <- getSeq( genome, foreground_granges_obj )
    names(foreground_sequences) <- foreground_granges_obj$SYMBOL#paste0("seq_", seq_along(foreground_sequences))
    writeXStringSet( foreground_sequences, "fg.fa")
    num_fg_sequences <- length(foreground_sequences)

    if (do_enrichment) {
        if (!is.null(background_granges_obj)) {
            background_sequences <- getSeq( genome, background_granges_obj )
            names(background_sequences) <- background_granges_obj$SYMBOL#paste0("seq_", seq_along(background_sequences))
        } else {
            background_sequences <- shuffle_sequences( foreground_sequences, method="euler", k=2)
        }
        writeXStringSet( background_sequences, "bg.fa")
        num_bg_sequences <- length(background_sequences)
    }

    tf_names <- names(pwm_list)
    num_motifs <- length(tf_names)
#    if (!enhancer_search_mode) {
#        motif_scores_mat <- matrix(0,nrow=num_fg_sequences,ncol=num_motifs,
#                              dimnames=list(names(foreground_sequences),tf_names))
#    } else {
        fg_motif_matches_mat <- matrix(0,nrow=num_fg_sequences,ncol=num_motifs,
                              dimnames=list(names(foreground_sequences),tf_names))
        bg_motif_matches_mat <- matrix(0,nrow=num_bg_sequences,ncol=num_motifs,
                              dimnames=list(names(background_sequences),tf_names))
#    }
    if (do_enrichment || enhancer_search_mode )
         enrichment_dt <- data.table( tf=tf_names )

    row_idx <- 1
    for (tf_name in tf_names) {
        print(paste(row_idx,length(tf_names),sep="/"))
        flush.console()
        pwm_obj <- pwm_list[[tf_name]]
        pwm_dt <- data.table( t(Matrix(merged_pwm_list[[tf_name]])) )
        pwm_file_name <- paste(tf_name,"pwm",sep=".")
        fwrite( data.table( pwm_dt ), pwm_file_name, col.names=F, row.names=F, quote=F, sep="\t" )

        #if (enhancer_search_mode) {
            fg_num_motif_matches <- get_sarus_motif_scores( "fg.fa", pwm_file_name,
            motif_p_value_thresh=motif_p_value_thresh )
            bg_num_motif_matches <- get_sarus_motif_scores( "bg.fa", pwm_file_name,
            motif_p_value_thresh=motif_p_value_thresh )

            if (!is.null(fg_num_motif_matches)) {
                num_fg_matches <- sum(fg_num_motif_matches)
                motif_containing_fg_seqs <- names(fg_num_motif_matches)
                fg_motif_matches_mat[motif_containing_fg_seqs,tf_name] <- fg_num_motif_matches
            } else {
                num_fg_matches <- 0
            }

            if (!is.null(bg_num_motif_matches)) {
                num_bg_matches <- sum(bg_num_motif_matches)
                motif_containing_bg_seqs <- names(bg_num_motif_matches)
                bg_motif_matches_mat[motif_containing_bg_seqs,tf_name] <- bg_num_motif_matches
            } else {
                num_bg_matches <- 0
            }
        #    enrichment_dt[row_idx,`:=`(num_fg=num_fg_matches,num_bg=num_bg_matches,fold_enrichment=enrichment)]
        #} else {
        #    fg_motif_scores <- get_sarus_motif_scores( "fg.fa", pwm_file_name, motif_p_value_thresh=motif_p_value_thresh )
        #    if (!is.null(fg_motif_scores)) {
        #        motif_scores_mat[names(fg_motif_scores),tf_name] <- fg_motif_scores
        #    }
        #}

        #if (do_enrichment && !enhancer_search_mode) {
        if (do_enrichment) {
            #bg_motif_scores <- get_sarus_motif_scores( "bg.fa", pwm_file_name )
            #motif_score_thresh <- quantile( bg_motif_scores, 1-bg_frac )
            #num_bg_matches <- floor( bg_frac * num_bg_sequences )

            #num_fg_matches <- sum(fg_motif_scores > motif_score_thresh)
            fg_frac_matched <- num_fg_matches/num_fg_sequences
            bg_frac_matched <- num_bg_matches/num_bg_sequences
            fisher_mat <- matrix(c(num_fg_matches,num_fg_sequences-num_fg_matches,
                                   num_bg_matches,num_bg_sequences-num_bg_matches),nrow=2,ncol=2,byrow=T)
            enr_p_value <- fisher.test( fisher_mat, alternative="greater" )$p.value

            #fold_enr <- fg_frac_matched/bg_frac
            fold_enr <- fg_frac_matched/bg_frac_matched

            enrichment_dt[row_idx,`:=`(fg_frac=fg_frac_matched,p_value=enr_p_value,fold_enrichment=fold_enr, bg_frac=bg_frac_matched,num_fg=num_fg_matches,num_bg=num_bg_matches)]
        }

        file.remove(pwm_file_name)

        row_idx <- row_idx + 1
    }
    return_list <- list()

#    if (!enhancer_search_mode) {
#        motif_scores_dt <- data.table( motif_scores_mat, keep.rownames=T ) %>% setnames(.,"rn","seq_name")
#        return_list[["features_dt"]] <- motif_scores_dt
#    } else {
#        print("Here")
#        4 + ""
        fg_motif_matches_dt <- data.table( fg_motif_matches_mat, keep.rownames=T ) %>% setnames(.,"rn","seq_name")
        bg_motif_matches_dt <- data.table( bg_motif_matches_mat, keep.rownames=T ) %>% setnames(.,"rn","seq_name")
        return_list[["fg_matches_dt"]] <- fg_motif_matches_dt
        return_list[["bg_matches_dt"]] <- bg_motif_matches_dt
#    }

    #if (do_enrichment || enhancer_search_mode) {
        #if (!enhancer_search_mode) {
            enr_q_values <- p.adjust( enrichment_dt$p_value, method="fdr" )
            enrichment_dt[,q_value:=enr_q_values]
        #}

        return_list[["enrichment_dt"]] <- enrichment_dt
    #}

    return(return_list)
}
