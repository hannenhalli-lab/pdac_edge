library(dplyr)
library(data.table)
library(ggplot2)
library(robustbase)
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

plot_skewness <- function(control_dist_dt) {
    combined_skewness_dt <- data.table()
    combined_skewness_dt <- data.table()
#    options(repr.plot.width=20, repr.plot.height=10)
#    theme_set(theme_gray(base_size = 20))
    for (cell_type_ in unique(control_dist_dt$cell_type)) {

        sub_dt <- control_dist_dt[cell_type == cell_type_]
        skewness_vals <- apply( sub_dt[,!c("cell_type","actual","cell.name")], 2, mc)
        skewness_dt <- data.table( sample_name=names(skewness_vals), skewness=skewness_vals )
        skewness_dt <- rbind( skewness_dt, data.table( sample_name="actual", skewness=mc(sub_dt$actual) ))

        control_mask <- grepl("control",skewness_dt$sample_name)
        skewness_dt$column_type <- "actual"
        skewness_dt[control_mask,]$column_type <- "control"
        actual_skewness <- skewness_dt[column_type == "actual",skewness]

        norm_fit_info <- gauss_fit(skewness_vals,actual_skewness,p_value_type="right.tail")

        #norm_p_value <- 1 - pnorm(actual_skewness,m,s )
        skewness_dt$cell_type <- paste( cell_type_, "p-value :", norm_fit_info$p_value )
        skewness_dt$norm_fit <- -1
        skewness_dt[control_mask,]$norm_fit <- norm_fit_info$norm_density
        combined_skewness_dt <- rbind( combined_skewness_dt, skewness_dt )
    }

    #colour_vec <- c("control"="black","actual"="red")
    p <- ggplot(data = combined_skewness_dt[column_type == "control",]) + geom_line( aes(x=skewness),stat="density",color="black") + geom_line(aes(x=skewness,y=norm_fit),linetype="dashed") + geom_vline( data = combined_skewness_dt[column_type == "actual",], aes(xintercept=skewness), color="red" ) +
    facet_wrap(~cell_type,ncol=3) + xlab("Skewness") + ylab("Density")

    return(p)
}

plot_edge_distance_ratio <- function(edge_malignant_dist_dt) {
    edge_malignant_dist_dt$plot_title <- ""
    edge_malignant_dist_dt$norm_fit <- -1
    for (cell_type_ in unique(edge_malignant_dist_dt$normal_cell_type) ) {
        sub_dt <- edge_malignant_dist_dt[normal_cell_type == cell_type_]
        num_controls <- nrow(sub_dt[dist_type == "control",])
        control_dists <- sub_dt[dist_type == "control",edge_center_dist_ratio]
        actual_dist <- sub_dt[dist_type == "actual",edge_center_dist_ratio] 
        norm_fit_info <- gauss_fit(control_dists,actual_dist, p_value_type="left.tail")

        plot_title <- paste( cell_type_, "p-value : ", norm_fit_info$p_value )
        edge_malignant_dist_dt[normal_cell_type == cell_type_,]$plot_title <- plot_title
        edge_malignant_dist_dt[normal_cell_type == cell_type_ & dist_type == "control",]$norm_fit <- norm_fit_info$norm_density
    }

    p <- ggplot( edge_malignant_dist_dt[dist_type=="control",] ) + geom_line(aes(x=edge_center_dist_ratio),stat="density",color="black") +
        geom_vline(data=edge_malignant_dist_dt[dist_type=="actual",],aes(xintercept=edge_center_dist_ratio),color="red") + geom_line(aes(x=edge_center_dist_ratio,y=norm_fit),color="black",linetype="dashed") + xlab("Edge distance ratio") + ylab("Density") + facet_wrap(~plot_title,ncol=3) + theme_classic()

    return(p)
}
