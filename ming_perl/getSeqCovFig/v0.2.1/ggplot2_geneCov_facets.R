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
source("geneCovPlot_fun_facets.R")
#covdir <- "~/Identify_sRNA/demo_H37Ra/0.Current_version/3.sRNA_viewer/covdir"
#gene <- "~/Identify_sRNA/demo_H37Ra/0.Current_version/3.sRNA_viewer/input/test.txt"
#image_fmt <- "png"
#outdir <- "test"
#gene <- "input/H37Rv_sRNA_list_100k_10k_1k_100.txt"

covidx <- getcovidx(covdir = covdir, chr = "H37Rv", readbig = FALSE)
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
      figs <- pos2fig(pos = pos, covidx = covidx, extend = 200) # defalut; gap = 200
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '.png')
        png(png_file, width = 1200, height = 1000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '.tiff')
        tiff(tiff_file, width = 1200, height = 1000, res = 300, compression = "zip+p")
      }
      print(figs[[1]])
      dev.off()
      ### free-y
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '-freey.png')
        png(png_file, width = 1200, height = 1000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '-freey.tiff')
        tiff(tiff_file, width = 1200, height = 1000, res = 300, compression = "zip+p")
      }
      print(figs[[2]])
      dev.off()
    }
    stopCluster(clu)
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}else if(image_fmt == "pdf"){
  ptime <- system.time({
#    pdf_file = paste0(outdir, '/', "seq_cov_plots.pdf")
    pdf_name <- gsub(".txt", ".pdf", basename(gene), fixed = TRUE)
    pdf_file <- paste0(outdir, "/", pdf_name)
    pdf(pdf_file, width = 6, height = 4)
    for(i in 1:nrow(genes)){
      pos  <- genes[i, ]
      figs <- pos2fig(pos = pos, covidx = covidx, extend = 200) # defalut; gap = 200
      ### Arrange multiple figs in one page
      print(figs[[1]])
      print(figs[[2]])
    }
    dev.off()
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}
