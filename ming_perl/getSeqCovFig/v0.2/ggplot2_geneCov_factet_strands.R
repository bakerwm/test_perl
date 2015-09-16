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
source("geneCovPlot_facet_strands.R")
covidx <- getcovidx(covdir = covdir, chr = "H37Rv", readbig = TRUE)
#gene <- "input/test.txt"
#gene <- "input/H37Rv_sRNA_list_100k_10k_1k_100.txt"
genes  <- readpos(gene)


library(bigmemory)
library(plyr)
library(ggplot2)
library(grid) # need for arrow plot
library(data.table)

### my theme
mythemeup <- theme_bw() +
  theme(#plot.title   = element_blank(),
    axis.line    = element_line(size = 0.7, colour = "black"),
    axis.text    = element_text(size = 8,   colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10,  colour = "black", 
                                face = "bold", vjust = 2.2),
    plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(colour = "black",
                                    size = 0.5, 
                                    fill = "gray70")
    )
mythemedown <- theme_bw() +
  theme(plot.title   = element_blank(),
        axis.line    = element_line(size = 0.7, colour = "black"),
        axis.text    = element_text(size = 8,   colour = "black"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,  colour = "black", 
                                    face = "bold", vjust = 2.2),
        plot.margin = unit(c(-0.4, 1, 0.5, 1), "cm"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(colour = "black",
                                        size = 0.5, 
                                        fill = "gray70")
        )
theme_arrow <- theme_bw() +
  theme(plot.title   = element_blank(),
        text         = element_text(size = 10, colour = "black"),
        plot.margin  = unit(c(-0.2, 1.5, -0.2, 2.2), "cm"),
        axis.ticks   = element_blank(),
        axis.line    = element_blank(),
        axis.text    = element_blank(),
        axis.title   = element_blank(),
        #panel.border = element_blank(),
        panel.grid   = element_blank()
        )

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
    da  <- foreach (i = 1:nrow(genes), .packages = c("ggplot2", "plyr", "scales", "grid", "bigmemory")) %dopar% { # For multiple CPUs
      #options( digits = 1 )
      pos  <- genes[i, ]
      figs <- pos2fig(pos = pos, covidx = covidx) # defalut; gap = 200
      p1 <- figs[[1]] + mythemeup #+ scale_y_continuous(labels = sci_10)
      p2 <- figs[[2]] + mythemeup #+ scale_y_continuous(labels = sci_10)
      p3 <- figs[[3]] + mythemedown #+ scale_y_continuous(labels = sci_10)
      p4 <- figs[[4]] + mythemedown #+ scale_y_continuous(labels = sci_10)
      ### Arrange multiple figs in one page
      library(gridExtra)
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '.png')
        png(png_file, width = 2000, height = 2000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '.tiff')
        tiff(tiff_file, width = 2000, height = 2000, res = 300, compression = "zip+p")
      }
      grid.arrange(p1, p3, ncol = 1, heights = c(9, 8))
      if(image_fmt %in% c("png", "tiff")){
        dev.off()
      }
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '-freey.png')
        png(png_file, width = 2000, height = 2000, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '-freey.tiff')
        tiff(tiff_file, width = 2000, height = 2000, res = 300, compression = "zip+p")
      }
      grid.arrange(p2, p4, ncol = 1, heights = c(9, 8))
      if(image_fmt %in% c("png", "tiff")){
        dev.off()
      }
    }
    stopCluster(clu)
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}else if(image_fmt == "pdf"){
  ptime <- system.time({
    pdf_name <- gsub(".txt", ".pdf", basename(gene), fixed = TRUE)
    pdf_file <- paste0(outdir, '/', pdf_name)
    #pdf_file = paste0(dirname(outdir), '/', "seq_cov_plots.pdf")
    pdf(pdf_file, width = 8, height = 8)
    for(i in 1:nrow(genes)){
      pos  <- genes[i, ]
      figs <- pos2fig(pos = pos, covidx = covidx) # defalut; gap = 200
      p1 <- figs[[1]] + mythemeup #+ scale_y_continuous(labels = sci_10)
      p2 <- figs[[2]] + mythemeup #+ scale_y_continuous(labels = sci_10)
      p3 <- figs[[3]] + mythemedown #+ scale_y_continuous(labels = sci_10)
      p4 <- figs[[4]] + mythemedown #+ scale_y_continuous(labels = sci_10)
      ### Arrange multiple figs in one page
      library(gridExtra)
      grid.arrange(p1, p3, ncol = 1, heights = c(9, 8))
      grid.arrange(p2, p4, ncol = 1, heights = c(9, 8))
    }
    dev.off()
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}
