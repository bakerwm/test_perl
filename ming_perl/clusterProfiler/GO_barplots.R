### bar plot for 
#setwd("/home/wangming/work/R_work/GO_analysis/wk_dirs/temp/GOdb_H37Rv")
GO_barplot <- function(x, 
                          showCategory = 5, 
                          enrich       = FALSE, 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1) {
  library(ggplot2)
  library(reshape2)
  library(grid)
  mytheme <- theme_bw() + 
    theme(text         = element_text(size = 8, colour = "black"),
          plot.title   = element_text(size = 10,  colour = "black",
                                      face = "bold", vjust = 0.7),
          axis.line    = element_line(size = 0.5, colour = "black"),
          axis.text    = element_text(size = 8,   colour = "black"),
          axis.title.x = element_text(size = 10,  colour = "black", 
                                      face = "bold", hjust = 0.5),
          axis.title.y = element_blank(),
          plot.margin  = unit(c(1,1,0.5,1), "cm"),
          panel.grid.minor  = element_blank())
  ### order the terms by Count, subset data
  fig_name <- gsub("\\.\\w+", "", basename(x), perl = TRUE)
  da <- read.table(x, sep = "\t", header = TRUE, quote = "")
  if(enrich) {
    dn <- subset(da, pvalue < pvalueCutoff & qvalue <= qvalueCutoff)
    if(nrow(dn) < showCategory) {
      showCategory <- nrow(dn)
    }
    de <- dn[ order(dn$pvalue), ]
    en <- de[1:showCategory, ]
    ### trim the name of Description
    en$Description <- rapply(as.list(en$Description), function(x) substring(x, 1, 40))
    en$Description <- as.factor(en$Description)
    ### delete unused levels
    #en$Description <- droplevels(en$Description)
    ### reorder the levels by pvalue
    en$Description <- reorder(en$Description, en$pvalue)
    #
    p <- ggplot(en, aes(x = Description, y = Count, fill = pvalue)) + 
      geom_bar(stat = "identity") +
      coord_flip() +
      ggtitle(fig_name) +
      #geom_text(aes(lable = Count, size = 8), hjust = 1.1, colour = "white") +
      scale_fill_gradient(low = "red", high = "blue") +
      mytheme 
  }else{
    da <- da[ order(-da$Count), ]
    if(nrow(da) < showCategory) {
      showCategory <- nrow(da)
    }
    group <- da[1:showCategory, ]
    ### delete the unused levels
    group$Description <- droplevels(group$Description)
    ### reorder the levels by count
    group$Description <- reorder(group$Description, group$Count, order = FALSE)
    #
    p <- ggplot(group, aes(x = Description, y = Count, fill = Count)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() + 
      ggtitle(fig_name) + 
      geom_text(aes(label = Count, size = 8), hjust = 1.1, colour = "white") +
      scale_fill_gradient(low="grey30", high="grey10") +
      mytheme + theme(legend.position = "none")
  }
  return(p)
}

#x <- "output_H37Rv/H37Rv_GO_category_BP.xls"
#x <- "output_H37Rv/H37Rv_GO_enrichment_BP.xls"
###
#usage <- "Usage: GO_barplot.R in.txt TRUE/FALSE(enrich)"
#args  <- commandArgs(trailingOnly = TRUE)
#if(length(args) != 2) stop(usage)
#f      <- args[1]
#enrich <- args[2]
#pdf_file <- gsub("\\.\\w+$", "_new.pdf", basename(f), perl = TRUE)
#
#pdf(pdf_file, width = 6, height = 4.5)
#GO_barplot(f, showCategory = 10, enrich = enrich, pvalueCutoff = 0.5)
#dev.off()
