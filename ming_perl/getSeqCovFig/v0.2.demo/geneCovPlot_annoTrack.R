### Generate coverage and gene symbols plot
setwd("~/Identify_sRNA/demo_H37Rv/0.Current_version/3.sRNA_viewer/")

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
        plot.margin = unit(c(1, 1, 0.1, 1), "cm")
  )
mythemedown <- theme_bw() +
  theme(plot.title   = element_blank(),
        axis.line    = element_line(size = 0.7, colour = "black"),
        axis.text    = element_text(size = 8,   colour = "black"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,  colour = "black", 
                                    face = "bold", vjust = 2.2),
        plot.margin = unit(c(-0.4, 1, 0.5, 1), "cm")
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

### cov plots
source("geneCovPlot_fun.R")
covidx <- getcovidx(covdir = "covdir", chr = "H37Rv", readbig = FALSE)
gene <- "input/test.txt"
#gene <- "input/H37Rv_sRNA_list_100k_10k_1k_100.txt"
genes  <- readpos(gene)

### ----parallel
library(parallel)
library(foreach)
library(doParallel)

pdf("H37Rv_sRNAs.pdf", width = 8, height = 8)
cores <- 12
clu   <- makeCluster(cores)
registerDoParallel(clu, cores = cores)
getDoParName()
ptime <- system.time({    
  da  <- foreach (i = 1:nrow(genes), .packages = c("ggplot2", "plyr", "scales", "grid", "bigmemory")) %dopar% { # For multiple CPUs
    options( digits = 1 )
    pos  <- genes[i, ]
    figs <- pos2fig(pos = pos, covidx = covidx) # defalut; gap = 200
    p1 <- figs[[1]] + mythemeup #+ scale_y_continuous(labels = sci_10)
    p2 <- figs[[2]] + mythemeup #+ scale_y_continuous(labels = sci_10)
    p3 <- figs[[3]] + mythemedown #+ scale_y_continuous(labels = sci_10)
    p4 <- figs[[4]] + mythemedown #+ scale_y_continuous(labels = sci_10)
    #p5 <- figs[[5]] + theme_arrow #+ scale_y_continuous(labels = sci_10)
    ### Arrange multiple figs in one page
    library(gridExtra)
    #pdf("test.pdf", width = 8, height = 8)
    #png_file <- paste0('output/images/', pos$id, '.png')
    #png(png_file, width = 2000, height = 2000, res = 300)
    #tiff_file <- paste0('output/images/', pos$id, '.tiff')
    #tiff(tiff_file, width = 2000, height = 2000, res = 300, compression = "zip+p")
    grid.arrange(p1, p3, ncol = 1, heights = c(9, 8))
    grid.arrange(p2, p4, ncol = 1, heights = c(9, 8))
    #grid.arrange(p2, p5, p4, ncol = 1, heights = c(8, 1, 7))
    #dev.off()
  }
  stopCluster(clu)
})
dev.off()

print(paste0("Finish (create pngs): ", ptime[3]))







