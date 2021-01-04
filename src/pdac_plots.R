library(dplyr)
library(data.table)
library(ggplot2)
library(robustbase)
library(stringr)
library(docxtools)

gauss_fit <- function(values_to_fit,test_val,p_value_type=NULL) {
    m <- mean(values_to_fit)
    s <- sd(values_to_fit)
    norm_density <- dnorm(values_to_fit,m,s)
    if (p_value_type == "right.tail") {
        norm_p_value <- 1 - pnorm(test_val,m,s )
    } else if (p_value_type == "left.tail") {
        norm_p_value <- pnorm(test_val,m,s )
    }

    return(list("p_value"=norm_p_value,"norm_density"=norm_density))
}

compute_skewness <- function(control_dist_dt) {
    combined_skewness_dt <- data.table()
    for (cell_type_ in unique(control_dist_dt$cell_type)) {
        sub_dt <- control_dist_dt[cell_type == cell_type_]
        control_skewness_vals <- apply( sub_dt[,!c("cell_type","actual","cell.name")], 2, mc)
        skewness_dt <- data.table( sample_name=names(control_skewness_vals), skewness=control_skewness_vals )
        skewness_dt <- rbind( skewness_dt, data.table( sample_name="actual", skewness=mc(sub_dt$actual) ))
        skewness_dt[,cell_type:=cell_type_]

        skewness_dt[,column_type:=ifelse(str_detect(sample_name,"control"),"control","actual")]
        actual_skewness <- skewness_dt[column_type == "actual",skewness]

        norm_fit_info <- gauss_fit(control_skewness_vals,actual_skewness,p_value_type="right.tail")

        skewness_dt[,p_value:=norm_fit_info$p_value]
        skewness_dt[,display_cell_type:=paste( cell_type_, "p-value :", norm_fit_info$p_value )]
        skewness_dt[column_type == "control",norm_fit:=norm_fit_info$norm_density]
        combined_skewness_dt <- rbind( combined_skewness_dt, skewness_dt )
    }

    return(combined_skewness_dt)
}

plot_skewness <- function(control_dist_dt,plot_facets=F,do_boxplot=F) {
    combined_skewness_dt <- compute_skewness(control_dist_dt)

    if (do_boxplot==T) {
            p <- ggplot(data = combined_skewness_dt[column_type == "control",]) + geom_line(aes(x=skewness),stat="density",linetype="dashed") + geom_vline( data = combined_skewness_dt[column_type == "actual",], aes(xintercept=skewness), size=1 ) + xlab("Skewness") + ylab("Density") + theme_pubr(base_size=15) + facet_wrap(~display_cell_type) + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
    } else {
        if (plot_facets) {
            p <- ggplot(data = combined_skewness_dt[column_type == "control",]) + geom_line(aes(x=skewness),stat="density",linetype="dashed") + geom_vline( data = combined_skewness_dt[column_type == "actual",], aes(xintercept=skewness), size=1 ) + xlab("Skewness") + ylab("Density") + theme_pubr(base_size=15) + facet_wrap(~display_cell_type) + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
        } else {
            p <- ggplot(data = combined_skewness_dt[column_type == "control",]) + geom_line(aes(x=skewness,color=cell_type),stat="density",linetype="dashed") + geom_vline( data = combined_skewness_dt[column_type == "actual",], aes(xintercept=skewness, color=cell_type), size=1 ) + xlab("Skewness") + ylab("Density") + theme_pubr(base_size=15) + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
        }
    }

    return(list("plot_obj"=p,"skewness_dt"=combined_skewness_dt))
}

compute_edge_distance_significance <- function(edge_malignant_dist_dt) {
    edge_malignant_dist_dt$plot_title <- ""
    edge_malignant_dist_dt$norm_fit <- -1
    for (cell_type_ in unique(edge_malignant_dist_dt$normal_cell_type) ) {
        sub_dt <- edge_malignant_dist_dt[normal_cell_type == cell_type_]
        num_controls <- nrow(sub_dt[dist_type == "control",])
        control_dists <- sub_dt[dist_type == "control",edge_center_dist_ratio]
        actual_dist <- sub_dt[dist_type == "actual",edge_center_dist_ratio] 
        norm_fit_info <- gauss_fit(control_dists,actual_dist, p_value_type="left.tail")

        facet_title <- paste( cell_type_, "p-value : ", norm_fit_info$p_value )
        edge_malignant_dist_dt[normal_cell_type == cell_type_,`:=`(plot_title=facet_title,p_value=norm_fit_info$p_value)]
        edge_malignant_dist_dt[normal_cell_type == cell_type_ & dist_type == "control",norm_fit:=norm_fit_info$norm_density]
    }

    return(edge_malignant_dist_dt)
}

plot_edge_distance_ratio <- function(edge_malignant_dist_dt,plot_facets=F,do_boxplot=F) {
    edge_malignant_dist_dt <- compute_edge_distance_significance(edge_malignant_dist_dt) 

    if (plot_facets) {
        p <- ggplot( edge_malignant_dist_dt[dist_type=="control",] ) +
        geom_line(aes(x=edge_center_dist_ratio),stat="density",color="black",linetype="dashed") +
        geom_vline(data=edge_malignant_dist_dt[dist_type=="actual",],aes(xintercept=edge_center_dist_ratio),color="red") + 
        xlab("Edge distance ratio") + ylab("Density") + facet_wrap(~plot_title,ncol=3) + 
        theme_pubr(base_size=15) + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
    } else {
        p <- ggplot( edge_malignant_dist_dt[dist_type=="control",] ) +
        geom_line(aes(x=edge_center_dist_ratio,color=normal_cell_type),stat="density",linetype="dashed") + 
        geom_vline(data=edge_malignant_dist_dt[dist_type=="actual",],aes(xintercept=edge_center_dist_ratio,color=normal_cell_type)) + xlab("Edge distance ratio") + ylab("Density") + theme_pubr(base_size=15) +
        theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
    }

    return(list("plot_obj"=p,"edge_distance_dt"=edge_malignant_dist_dt))
}

plot_tumourwise_fingerprints <- function(tumour_wise_edgeness_info_list) {
    tumour_samples <- names(tumour_wise_edgeness_info_list)
    combined_skewness_dt <- data.table()
    combined_edge_distance_dt <- data.table()

    for (tumour_sample in tumour_samples) {
        edge_info_list <- tumour_wise_edgeness_info_list[[tumour_sample]]
        skewness_dt <- compute_skewness(edge_info_list$control_dist)
        skewness_dt <- skewness_dt[,tumour:=tumour_sample] %>% dplyr::select(tumour,cell_type,column_type,p_value,skewness) 
        combined_skewness_dt <- rbind(combined_skewness_dt,skewness_dt)
        edge_distance_dt <- compute_edge_distance_significance( edge_info_list$edge_malignant_dist_dt )
        edge_distance_dt[,tumour:=tumour_sample]
        combined_edge_distance_dt <- rbind(combined_edge_distance_dt,edge_distance_dt)
    }
    combined_skewness_dt <- combined_skewness_dt %>% group_by(cell_type,tumour,p_value) %>% summarize() %>% mutate(p_adj=p.adjust(p_value,method="bonferroni")) %>% dplyr::select(cell_type,tumour,p_adj) %>% merge(combined_skewness_dt,by=c("cell_type","tumour")) %>% mutate(significant=ifelse(p_adj < 0.1,TRUE,FALSE)) %>% data.table
    tumour_sample_order <- paste0("T",1:24)
    combined_skewness_dt[,tumour:=factor(tumour,levels=tumour_sample_order)]

    combined_edge_distance_dt <- combined_edge_distance_dt %>% group_by(normal_cell_type,tumour,p_value) %>% summarize() %>% mutate(p_adj=p.adjust(p_value,method="bonferroni")) %>% dplyr::select(normal_cell_type,tumour,p_adj) %>% merge(combined_edge_distance_dt,by=c("normal_cell_type","tumour")) %>% mutate(significant=ifelse(p_adj < 0.1,TRUE,FALSE)) %>% data.table
    combined_edge_distance_dt[,tumour:=factor(tumour,levels=tumour_sample_order)]

    skewness_plot_obj <- ggplot( combined_skewness_dt[column_type == "control"] ) +
    geom_boxplot(aes(x=tumour,y=skewness),color="gray",fill=NA) + geom_point(data=combined_skewness_dt[column_type ==
    "actual",],aes(x=tumour,y=skewness,color=significant),size=4) + theme_pubr(base_size=20) + theme(axis.text.x=element_text(angle=90,vjust=0.5)) + facet_wrap(~cell_type,ncol=1) + ylab("Skewness") + xlab(NULL) + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))
    print(skewness_plot_obj)

    edge_distance_plot_obj <- ggplot( combined_edge_distance_dt[dist_type == "control"] ) +
    geom_boxplot(aes(x=tumour,y=edge_center_dist_ratio),color="gray",fill=NA) +
    geom_point(data=combined_edge_distance_dt[dist_type ==
    "actual",],aes(x=tumour,y=edge_center_dist_ratio,color=significant),size=4) + theme_pubr(base_size=20) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5)) + facet_wrap(~normal_cell_type,ncol=1) + xlab(NULL) + ylab("Edge/Non-Edge
    Distance Ratio")  + theme(panel.grid.major=element_line(color="gray",linetype="dotted"))

    print(edge_distance_plot_obj)
}
