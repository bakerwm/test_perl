### ----This file contains three functions for GO analysis-----------------
### 1.getGOmap: buildGOmap() for your GO annotation files
### 2.gene2GO: for default GO analysis in clusterProfiler
### 3.GO_barplot: customized barplots for GO analysis
###
### Wang Ming wangmcas(AT)gmail.com
### 2015-07-20

### perform GO analysis using clusterProfiler package
getGOmap <- function(godb) {
  suppressPackageStartupMessages(library(clusterProfiler))
  ### build GOmap if required
  ### two column file, GO annotation (column 1) and Entrez GeneID (column 2)
  #godb <- "/home/wangming/work/R_work/GO_analysis/wk_dirs/temp/GOdb_H37Rv/H37Rv.GO.db"
  ### build
  gomap <- read.table(godb, colClasses = "character",
                      col.names = c("go_accession", "entrezgene"))
  ### check input data (begin with GO annotation)
  if( nrow(gomap) != sum(grepl('GO:', gomap[, 1]))) {
    f1 <- "GO.db file is not in this format: <GO accession> <GeneID>"
    f  <- paste0(godb, w2, sep = " ")
    stop(f)
  }
  dim(gomap)
  head(gomap, 4)
  buildGOmap(gomap)
}

gene2GO <- function(gene, gff, godb, gomap = FALSE, strain = "H37Rv", outdir = "output", ...) {
  # gene, character
  # gff,  file NCBI format
  # godb, <GO> <GeneID>
  # strain, the title
  # outdir, the dir for GO analysis ouput
  ### load package
  suppressPackageStartupMessages(library(DOSE))
  suppressPackageStartupMessages(library(clusterProfiler))
  # create output dir
  #outdir <- paste("output", strain, sep = "_")
  if(! outdir %in% list.files("./") ) {
    dir.create(outdir, showWarnings = FALSE, mode = "0755")
  }
  ### load GFF file, to transform the locus tag to GeneID
  #gff <- "/home/wangming/work/database/H37Rv/NC_000962.gff"
  Gff2GeneTable(gff)
  load("geneTable.rda")
  head(geneTable, 4)
  # perform GO analysis
  ### buildGOmap ?, Y/N
  if(gomap) {
    getGOmap(godb)
  }
  # prepare data
  gene_id <- as.character(geneTable[geneTable$Locus %in% gene, "GeneID"])
  ### prepare input data
  if( ! is.character(gene_id) ) {
    stop("Input gene need to be: character")
  }
  ### Start GO analysis
  ### 1. GOgroup(), classification
  ### organism is not specified, have to not included in default 20 organisms
  ont  <- c("BP", "CC", "MF")
  type <- c("_GO_category_", "_GO_enrichment_", "_GO_gse_")
  fmt  <- c(".xls", ".pdf")
  ont  <- as.list(ont)
  fpre <- rapply(ont, function(x) paste0(strain, type, x))
  xls  <- paste(outdir, "/", fpre, ".xls", sep = "")
  fig  <- paste(outdir, "/", fpre, ".pdf", sep = "")
  for(i in 1:length(ont)) {
    ### 1. GO categories
    ggo1   <- groupGO(gene_id, organism = strain, ont = ont[[i]], level = 3, readable = TRUE)
    da1    <- summary(ggo1)
    da1    <- da1[ order(-da1[, 3]), ]
    head(da1, 4)
    bp_ggo <- barplot(ggo1, drop = TRUE, showCategory = 10)
    write.table(da1, xls[3*i - 2], append = FALSE, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    pdf(fig[3*i - 2], width = 8, height = 6)
    print(bp_ggo)
    dev.off()
    ### 2. GOenrich(), enrichment analysis
    ego1   <- enrichGO(gene_id, organism = strain, ont = ont[[i]], 
                      pvalueCutoff = 1, qvalue = 1, readable = TRUE)
    da2    <- summary(ego1)
    da2    <- da2[ order(-da2[, "Count"]), ]
    bp_ego <- barplot(ego1, drop = TRUE, showCategory = 10)
    write.table(da2, xls[3*i - 1], append = FALSE, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    pdf(fig[3*i - 1], width = 20, height = 16)
    print(bp_ego)
    enrichMap(ego1)
    cnetplot(ego1, categorySize = "pvalue")
    dev.off()
    ### 3. GO Gene Set Enrichment Analysis
    #ego2   <- gseGO(gene_id, organism = strain, ont = ont[[i]], nPerm = 1, minGSSize = 1,
    #              pvalueCutoff = 1, verbose = FALSE)
    #da3    <- summary(ego2)
    #write.table(da3, xls[3*i - 0], append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    #pdf(fig[3*i - 0], width = 10, height = 8)
    #barplot(ego2, drop = TRUE, showCategory = 10)
    #enrichMap(ego2)
    #cnetplot(ego2, categorySize = "pvalue")
    #dev.off()
  } 
  ### 4. KEGG analysis
  ### Find the support organism in the full list:
  ### http://www.genome.jp/kegg/catalog/org_list.html
  ### Common use: 
  ### mtu H37Rv (1998)
  ### mtv H37Rv (2012)
  ### mra H37Ra (2007)
  ### mtc CDC1551 (2001)
  ### msm Msmeg_MC2_155 (2006)
  ### msg Msmeg_MC2_155 (2012)
  ### msg Msmeg_MC2_155 (2014)
  ### mbb Mbovis_BCG_Pasteur (2007)
  ### convert GeneID to locus name
#  gene_locus <- geneTable[geneTable$GeneID %in% gene_id, "Locus"]
#  ### (1) KEGG over-representation test
#  kk1 <- enrichKEGG(gene_locus, organism = "mtu", pvalueCutoff = 1, qvalueCutoff = 1)
#  dk1 <- summary(kk1)
#  ### (2) KEGG Gene Set Enrichment Analysis
#  #kk2 <- gseKEGG(gene_locus, organism = "mtu", nPerm = 1, minGSSize = 2, pvalueCutoff = 1, verbose = FALSE)
#  #head(summary(kk2))
#  kk_xls <- paste0(outdir, "/", strain, "_KEGG_enrichment.xls")
#  kk_pdf <- paste0(outdir, "/", strain, "_KEGG_enrichment.pdf")
#  write.table(dk1, kk_xls, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#  pdf(kk_pdf, width = 20, height = 16)
#  enrichMap(kk1)
#  cnetplot(kk1, categorySize = "geneNum")
#  dev.off()
#  xlslist <- c(xls, kk_xls)
   xlslist <- xls
   return(xlslist)
}

### ----Recreate the barplots, custemised theme()--------------------------
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

### ----END OF FILE--------------------------------------------------------
