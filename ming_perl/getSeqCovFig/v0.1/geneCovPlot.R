### ----plot the sRNA on reference genome----------------------------------
### ----four libraries
### ----parse coverage file, index-----------------------------------------
getcovidx <- function(covdir, chr = "H37Rv", readbig = FALSE) {
  covfiles <- list.files(path = covdir, pattern = ".coverage.\\w$", full.names = TRUE)
  if(length(covfiles) < 2) stop(paste0("cov files not found in: ", covdir))
  df <- data.frame(do.call(rbind, strsplit(basename(covfiles), split = '.', fixed = TRUE)))
  ### lib type is 2-digit type
  df$libtype <- substr(df$X1, nchar(chr) + 1, nchar(chr) + 2)
  covidx <- data.frame(chr <- chr, libtype <- df$libtype, strand <- df$X3)
  names(covidx) <- c("chr", "libtype", "strand")
  covidx$file.path <- covfiles
  ### local *.des and *.bin files
  covidx$bin <- paste0(covfiles, '.bin')
  covidx$desc <- paste0(covfiles, '.desc')
  ###
  if(readbig) {
    tp <- mapply(readbigfile, covidx$file.path, covidx$bin, covidx$desc, covdir)
  }
  ### convert the strand name to: +, -
  #library(dplyr)
  lut <- c("n" = "-", "p" = "+")
  covidx$strand <- lut[covidx$strand]
  return(covidx)
}
###
readbigfile <- function(x, bin, desc, covdir) {
  library(bigmemory)
  read.big.matrix(x, sep = "\t", header = FALSE, type = "integer",
                  col.names = c("chr", "position", "count"),
                  backingpath = dirname(bin),
                  backingfile = basename(bin),
                  descriptor  = basename(desc))
}
###
pos2cov <- function(x, covidx, gap = 200) {
  subidx <- subset(covidx, strand == x$strand)
  subcov <- data.frame()
  for(i in 1:nrow(subidx)) {
    ### define the range of the x-axis:
    start    <- x$start - gap
    if(start < 0) start <- 1
    end      <- x$end + gap
    ###
    df1      <- attach.big.matrix(subidx$desc[i], subidx$bin[i]) # desc, bin
    df2      <- df1[mwhich(df1, c("position"), list(c(start, end)),
                          list(c('ge', 'le')), 'AND')]
    df2      <- as.data.frame(df2)
    df2$id   <- x$id
    df2$chr  <- x$chr
    df2$libtype <- subidx$libtype[i]
    rm(start, end, df1)
    subcov <- rbind.fill(subcov, df2)
  }
  return(subcov)
}
###
pos2fig   <- function(pos, covidx, gap = 200) {
  #pos
  subcov  <- pos2cov(pos, covidx, gap = 200)
  ### axis
  xstart  <- ceiling(min(subcov$position) / 10) * 10
  xend    <- floor(max(subcov$position) / 10) * 10
  xstep   <- floor(floor((xend - xstart) / 7) / 10) * 10
  xbreaks <- seq(xstart, xend, by = xstep)
  ###
  ### change lib name
  libname <- c("01" = "18-40 nt", "02" = "40-80 nt",
               "03" = "80-140 nt", "04" = ">140 nt")
  subcov$libtype <- as.factor(libname[subcov$libtype])
  subcov$libtype <- factor(subcov$libtype, levels = c("18-40 nt", "40-80 nt", "80-140 nt", ">140 nt"))
  ###
  mytheme <- theme_bw() +
    theme(plot.title = element_text(colour = "black", size = 20),
          axis.line  = element_line(colour = "black", size = .8),
          axis.title = element_text(colour = "black", size = 16),
          axis.text  = element_text(colour = "black", size = 12),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(colour = "black",
                                          size = 0.5, 
                                          fill = "gray70"))
  p <- ggplot(subcov, aes(x = position, y = count)) +
    geom_area(fill = "grey20") +
    geom_vline(xintercept = c(pos$start, pos$end), colour = "blue",
               size = .2, linetype = "longdash") +
    scale_x_continuous(breaks = xbreaks) +
    ggtitle(pos$note) + 
    mytheme
  p1 <- p + facet_grid(libtype ~ .)
  p2 <- p + facet_grid(libtype ~ ., scales = "free_y")
  return(list(p1, p2))
}

readpos <- function(x) {
  df <- read.table(x, header = FALSE, sep = "\t", na.strings = "NA")
  df2 <- df[, 1:6]
  names(df2) <- c('id', 'chr', 'length', 'start', 'end', 'strand')
  df2 <- within(df2,  note <- paste0(chr, '_', id, ' (', length, ' nt) ', start, ':', end, ':', strand))
}

### ----main script--------------------------------------------------------
#setwd("/home/wangming/Identify_sRNA/demo_H37Rv/0.Current_version/3.sRNA_viewer")
#################################
# Load parameters
Usage <- "Usage: Rscript getSeqCovFig.R <covdir> <gene> <outdir> <readcov|TRUE,FALSE>"
args  <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop(Usage)

covdir <- args[1]
gene   <- args[2]
outdir <- args[3]
readcov <- args[4]

pngdir <- paste(outdir, "images", sep = '/')
#if(dir.exists(outdir)) stop(paste0("Need a new folder: ", outdir))
if(! dir.exists(pngdir)) dir.create(pngdir, showWarnings = TRUE, recursive = TRUE)
print(paste(c("Input para:", args), collapse=" "))

### ----Read input file----------------------------------------------------
covidx <- getcovidx(covdir = "covdir", chr = "H37Rv", readbig = readcov)
genes  <- readpos(gene)

### ----parallel
library(parallel)
library(foreach)
library(doParallel)

cores <- 10
clu   <- makeCluster(cores)
registerDoParallel(clu, cores = cores)
getDoParName()
ptime <- system.time({    
  da  <- foreach (i = 1:nrow(genes), .packages = c("ggplot2", "plyr", "bigmemory")) %dopar% { # For multiple CPUs
    pos  <- genes[i, ]
    figs <- pos2fig(pos = pos, covidx = covidx) # defalut; gap = 200
    pngfile <- paste0(pngdir, "/", as.character(pos$id), c(".png", "-freey.png"))
    #tifffile <- paste0(pngdir, '/', as.character(pos$id), c(".tiff", "-freey.tiff"))
    png(pngfile[1], width = 800, height = 600)
    #tiff(pngfile[1], width = 7, height = 5, units = "in", res = 150)
    print(figs[[1]])
    dev.off()
    png(pngfile[2], width = 800, height = 600)
    #tiff(pngfile[2], width = 7, height = 5, units = "in", res = 150)
    print(figs[[2]])
    dev.off()
  }
  stopCluster(clu)
})
print(paste0("Finish (create pngs): ", ptime[3]))

### generate htmls for png files
#cmd1 <- paste("mv", OutDir, "images", sep = " ")
#cmd2 <- paste("mv", "images", OutDir, sep = " ")
#cmd3 <- paste("perl", "png_to_html.pl", InFile, OutDir, sep = " ")
#
#t1 <- try(system(cmd1, intern = TRUE))
#dir.create(OutDir, showWarnings = TRUE)
#t2 <- try(system(cmd2, intern = TRUE))
#t3 <- try(system(cmd3, intern = TRUE))
