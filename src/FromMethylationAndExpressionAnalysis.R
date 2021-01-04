library(ggplot2)
library(ggpubr)
ggTest <- list( c("EdgeUp", "RestOfGenome"), c("EdgeDown", "RestOfGenome"))



########################Expression@##########################
medianTumorZscores <- read.table("MedianZscoreOfTumorSamplesForEdgeCellGenes.txt", sep = "\t", header = T)

medianTumorZscores$GeneType <- factor(medianTumorZscores$GeneType, levels = c("EdgeUp", "RestOfGenome", "EdgeDown"))
maxLimit <- boxplot.stats(medianTumorZscores$Zscore)$stats[c(1, 5)]	

bpG <- ggboxplot(na.omit(medianTumorZscores), 
				x="GeneType", y="Zscore",
				ylab=("Tumor Z-score"),
				color = "GeneType", lwd = 0.5, outlier.shape = NA, facet.by = "CellType")
				
bpG <- bpG  + 	theme_bw() +
			theme(axis.text.x = element_blank(),
				axis.text.y = element_text(size = 10),
				axis.title.x = element_blank(),
				legend.position = "bottom",
				axis.title.y = element_text(size = 10),
				axis.ticks.x = element_blank())

ggsave( "TumorZscoresForEdgeCellGenesFromGtex.svg", bpG + coord_cartesian(ylim = (maxLimit + c(-30, 7))) + 
	stat_compare_means(comparisons = ggTest, label = "p.format", size = 3, label.y = c(10,11)),
width = 5, height = 3, pointsize = 10, units = 'in' )


##############################Methylation################################
methylationDf <- read.table("TumorAndNoralSamplesMethylationForEdgeCellGenes.txt", sep = "\t", header = T)


methylationDf$GeneType <- factor(methylationDf$GeneType, levels = c("EdgeUp", "RestOfGenome", "EdgeDown"))
#library(dplyr)

#############Plot delta##############
maxLimit <- boxplot.stats(methylationDf$deltaMethylation)$stats[c(1, 5)]	
bpG <- ggboxplot(na.omit(methylationDf), 
				x="GeneType", y="deltaMethylation",
				ylab=("Delta (Tumor - Normal)"),
				color = "GeneType", lwd = 0.5, outlier.shape = NA, facet.by = "CellType")
				
bpG <- bpG  + 	theme_bw() +
			theme(axis.text.x = element_blank(),
				#axis.text.x = element_text(size = 12, face = "bold"),
				axis.text.y = element_text(size = 10),
				axis.title.x = element_blank(),
				legend.position = "bottom",
				axis.title.y = element_text(size = 10))


ggsave("MethylationDeltaAmongGeneCategoriesFromArrayExpress.svg",
bpG + coord_cartesian(ylim = (maxLimit + c(-0.2, 0.3))) + 
	stat_compare_means(method = "wilcox.test", comparisons = ggTest, label = "p.format", size = 3, label.y =
    c(0.3,0.35)), width = 5, height = 3, pointsize = 10, units = 'in')
