
### ----ClusterPorfiler default analysis-----------------------------------
### load functions
setwd("/home/wangming/work/R_work/GO_analysis/wk_dirs/Rv")
source("gene2GOanalysis.R")
gff  <- "/home/wangming/work/database/H37Rv/NC_000962.gff"
godb <- "/home/wangming/work/R_work/GO_analysis/QuickGO_annotation/GO_annotation/NC_000962.GO.db"

### for down
load("20150723_sRNA_targets_down.rda")
gene <- gsub("rv", "Rv", tgdown$locus, ignore.case = FALSE)
filelist <- gene2GO(gene = gene, gff = gff, godb = godb, outdir = "tar_down", gomap = FALSE)

### for up
load("20150723_sRNA_targets_up.rda")
gene <- gsub("rv", "Rv", tgup$locus, ignore.case = FALSE)
filelist <- gene2GO(gene = gene, gff = gff, godb = godb, outdir = "tar_up", gomap = FALSE)

### ----Recreate barplots()------------------------------------------------
calist <- filelist[grepl("category", filelist, perl = TRUE)]
enlist <- filelist[grepl("enrichment", filelist, perl = TRUE)]
make_plot <- function(x, enrich = FALSE, showCategory = 10, pvalueCutoff = 0.5) {
  pdf_out <- gsub("\\.\\w+$", "_new.pdf", x, perl = TRUE)
  pdf(pdf_out, width = 5, height = 4.5)
  p <- GO_barplot(x, showCategory = showCategory, enrich = enrich, 
                  pvalueCutoff = pvalueCutoff)
  print(p)
  dev.off()
}
sapply(calist, function(x) make_plot(x, showCategory = 10, enrich = FALSE))
sapply(enlist, function(x) make_plot(x, showCategory = 10, enrich = TRUE, pvalueCutoff = 0.5))

### ----H37Rv H37Ra mRNA differential Expression---------------------------
mRNA_diff <- "/home/wangming/work/Rv_Ra_sRNA/exp_analysis/mRNA/H37Ra_H37Rv_exp.diff"
dd        <- read.table(mRNA_diff, header = TRUE, sep = "\t")
dd2       <- dd[c("test_id", "value_1", "value_2", "log2.fold_change.")]
mRdown    <- dd2$test_id[dd2$log2.fold_change. < -1]
mRup      <- dd2$test_id[dd2$log2.fold_change. > 1]

gene <- as.character(mRdown)
filelist <- gene2GO(gene = gene, gff = gff, godb = godb, outdir = "mR_down", gomap = FALSE)

gene <- as.character(mRup)
filelist <- gene2GO(gene = gene, gff = gff, godb = godb, outdir = "mR_up", gomap = FALSE)

### ----END OF FILE--------------------------------------------------------
