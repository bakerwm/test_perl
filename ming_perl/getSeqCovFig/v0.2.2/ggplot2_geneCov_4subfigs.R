### Generate coverage and gene symbols plot
#setwd("~/Identify_sRNA/demo_H37Rv/0.Current_version/3.sRNA_viewer/")

#################################
# Load parameters
Usage <- "Usage: Rscript ggplot2_geneCov.R <covDir> <gene.txt> <out.fmt/png,tiff,pdf> <outdir>"
args  <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop(Usage)
covdir <- args[1]
gene   <- args[2]
image_fmt <- args[3]
outdir <- args[4]
print( paste(c("Input para:", args), collapse=" ") )
if(! dir.exists(outdir)){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}
#################################

### cov plots
source("geneCovPlot_fun_4subfigs.R")
covidx <- getcovidx(covdir = covdir, chr = "H37Rv", readbig = FALSE)
#gene <- "input/test.txt"
#gene <- "input/H37Rv_sRNA_list_100k_10k_1k_100.txt"
genes  <- readpos(gene)


library(bigmemory)
library(plyr)
library(ggplot2)
library(grid) # need for arrow plot
library(data.table)

### my theme

### ----parallel-----------------------------------------------------------
if(image_fmt %in% c("png", "tiff")){
  library(parallel)
  library(foreach)
  library(doParallel)
  
  cores <- 16
  clu   <- makeCluster(cores)
  registerDoParallel(clu, cores = cores)
  getDoParName()
  ptime <- system.time({    
    da  <- foreach (i = 1:nrow(genes), .packages = c("ggplot2", "plyr", "scales", "grid", "gridExtra", "bigmemory")) %dopar% { # For multiple CPUs
      #options( digits = 1 )
      pos  <- genes[i, ]
      library(gridExtra)
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '.png')
        png(png_file, width = 2000, height = 1000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '.tiff')
        tiff(tiff_file, width = 2000, height = 1000, res = 300, compression = "zip+p")
      }
      figs <- pos2fig(pos = pos, covidx = covidx, yfree = FALSE) # defalut; gap = 200
      p1 <- figs[[1]] + theme(plot.margin = unit(c(0.8, 0.1, -0.1, 0.2), "cm"))
      p2 <- figs[[2]] + theme(plot.margin = unit(c(0.8, 0.5, -0.1, 0.1), "cm"))
      p3 <- figs[[3]] + theme(plot.margin = unit(c(-0.1, 0.1, 0.1, 0.2), "cm"))
      p4 <- figs[[4]] + theme(plot.margin = unit(c(-0.1, 0.5, 0.1, 0.1), "cm"))
      grid.arrange(p1, p2, p3, p4, ncol = 2,
                   top = pos$note, left = "Count", bottom = "Position")
      if(image_fmt %in% c("png", "tiff")){
        dev.off()
      }
      ### free-y
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '-freey.png')
        png(png_file, width = 2000, height = 1000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '-freey.tiff')
        tiff(tiff_file, width = 2000, height = 1000, res = 300, compression = "zip+p")
      }
      figs <- pos2fig(pos = pos, covidx = covidx, yfree = TRUE) # defalut; gap = 200
      p1 <- figs[[1]] + theme(plot.margin = unit(c(0.8, 0.1, -0.1, 0.2), "cm"))
      p2 <- figs[[2]] + theme(plot.margin = unit(c(0.8, 0.5, -0.1, 0.1), "cm"))
      p3 <- figs[[3]] + theme(plot.margin = unit(c(-0.1, 0.1, 0.1, 0.2), "cm"))
      p4 <- figs[[4]] + theme(plot.margin = unit(c(-0.1, 0.5, 0.1, 0.1), "cm"))
      grid.arrange(p1, p2, p3, p4, ncol = 2,
                   top = pos$note, left = "Count", bottom = "Position")
      if(image_fmt %in% c("png", "tiff")){
        dev.off()
      }
    }
    stopCluster(clu)
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}

if(image_fmt == "pdf"){
  ptime <- system.time({
    pdf_file = paste0(outdir, '/', "seq_cov_plots.pdf")
    pdf(pdf_file, width = 10, height = 5)
    for(i in 1:nrow(genes)){
      pos  <- genes[i, ]
      figs <- pos2fig(pos = pos, covidx = covidx, yfree = FALSE) # defalut; gap = 200
      p1 <- figs[[1]] + theme(plot.margin = unit(c(0.8, 0.1, -0.1, 0.2), "cm"))
      p2 <- figs[[2]] + theme(plot.margin = unit(c(0.8, 0.5, -0.1, 0.1), "cm"))
      p3 <- figs[[3]] + theme(plot.margin = unit(c(-0.1, 0.1, 0.1, 0.2), "cm"))
      p4 <- figs[[4]] + theme(plot.margin = unit(c(-0.1, 0.5, 0.1, 0.1), "cm"))
      ### Arrange multiple figs in one page
      library(gridExtra)
      grid.arrange(p1, p2, p3, p4, ncol = 2,
                   top = pos$note, left = "Count", bottom = "Position")
    }
    dev.off()
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}
