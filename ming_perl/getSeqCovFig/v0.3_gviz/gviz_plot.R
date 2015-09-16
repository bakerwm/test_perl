# Using Gviz package to generate RNA-seq reads coverage maps, combined with sequence symbols (Gene, sRNA)
#setwd("/home/wangming/work/R_work/gviz_plots")

# The Gviz plots consist of multiple tracks, usually, Axis track, 
# ideogram track, Data Track, Annotation Track, Data Track, and 
# gene track

#################################
# Load parameters
Usage <- "Usage: Rscript gviz_plot.R <covDir> <gff> <gene.bed> <out.fmt/png,tiff,pdf> <outdir>"
args  <- commandArgs(trailingOnly=TRUE)
if(length(args) != 6) stop(Usage)
covdir    <- args[1]
gff_file  <- args[2]
sR_file   <- args[3]
sR_file2  <- args[4]
image_fmt <- args[5]
outdir    <- args[6]
print( paste(c("Input para:", args), collapse=" ") )
if(! dir.exists(outdir)){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}
#################################

### ----prepare data-------------------------------------------------------
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
#library(GenomicFeatures)

## 1.GenomeAxisTrack
axisTrack <- GenomeAxisTrack(labelPos     = "alternating",
                             #distFromAxis = 0.5,
                             fontsize     = 9,
                             fontcolor    = "black",
                             showId       = TRUE,
                             col.id       = "black",
                             cex.id       = 0.7,
                             fill.range   = "lightgreen",
                             exponent     = 0)
#plotTracks(axisTrack, from = gfrom, to = gto)

## 2.1 Genome Annotation gene track
#gff_file <- "input/NC_000962.gff"
gff_gr   <- import.gff(gff_file)
gff_gr   <- gff_gr[gff_gr@elementMetadata$type == "gene", 
                   c("locus_tag", "type", "gene")]
names(values(gff_gr)) <- c("group", "type", "name")

## 2.2 sRNA Track
#sR_file  <- "input/H37Rv_sRNA.bed"
#sR_file  <- "input/test.bed"
sR_gr    <- import.bed(sR_file)  # name, score
names(values(sR_gr)) <- c("id", "length")
sR_df    <- as(sR_gr, "data.frame")

## 2.3 other seq Track (.bed)
sR_gr2   <- import.bed(sR_file2)
names(values(sR_gr2)) <- c("id", "length")
sR_df2   <- as(sR_gr2, "data.frame")

## 2.1 AnnotationTrack
options(ucscChromosomeNames=FALSE)
annTk <- function(range, name = "Genome", fill = "grey70", anno = "id") {
  atk <- AnnotationTrack(range, 
                         name              = name, 
                         fill              = fill,
                         shape             = "fixedArrow", 
                         stacking          = "full", 
                         groupAnnotation   = anno,
                         fontcolor.group   = "darkblue",
                         fontface.group    = 2,
                         fontsize.group    = 10,
                         just.group        = "below",
                         arrowHeadMaxWidth = 1,
                         min.width         = 3, 
                         min.height        = 10,
                         collapse          = TRUE, 
                         showOverplotting  = TRUE)
}
aTrack  <- annTk(gff_gr, name = "Genes", fill = "grey70", anno = "group")
ncTrack <- annTk(sR_gr, name = "sRNAs", fill = "deepskyblue", anno = "id")
sTrack  <- annTk(sR_gr2, name = "sRNAs", fill = "lightgreen", anno = "id")
#plotTracks(aTrack, from = gfrom, to = gto)
#plotTracks(ncTrack, from = gfrom, to = gto)

## 5.DataTrack
## bedGraph, BED, ...
##cov_file  <- "input/H37Rv01.fwd.bedgraph"
##cov_gr    <- import.bedGraph(cov_file)
cov_files <- list.files(covdir, "*.bedgraph", full.names = TRUE)
cov <- lapply(cov_files, function(x) import.bedGraph(x))

### ----Function: create DataTrack-----------------------------------------
dtk <- function(range, 
                name      = "Data",
                fill_hist = "darkgrey"){
  dtrack <- DataTrack(range, 
                      name           = name, 
                      type           = c("histogram", "g"),
                      fill.histogram = fill_hist,
                      lty            = "blank",
                      lwd            = 0.4,
                      lty.grid       = 0,
                      window         = -1,
                      windowSize     = 1,
                      cex.axis       = 0.6,
                      col.axis       = "grey20")
  return(dtrack)
}

### prepare library index
#dtk1p <- dtk(cov_1p, name = "Lib01 +", fill_hist = "red")
fill_hist <- rep(c("red", "blue"), times = 4)
name      <- paste0("Lib0", rep(1:4, each = 2), " ", rep(c("+", "-")))
dTk       <- lapply(1:8, function(x) dtk(cov[[x]], name = name[x], fill_hist = fill_hist[x]))

# combine all tracks
#pdf("demo1.pdf", width = 12, height = 6)
pts <- function(gfrom, gto, main, ymax = 1e3, extend = 200) {
  plotTracks(list(dTk[[7]], dTk[[5]], dTk[[3]], dTk[[1]], 
                  dTk[[8]], dTk[[6]], dTk[[4]], dTk[[2]], 
                  axisTrack, aTrack, ncTrack, sTrack),
             sizes            = c(rep(1, 8), 1, 1.1, 1.1, 1.1),
             from             = gfrom - extend,
             to               = gto + extend,
             main             = main,
             ylim             = c(0, ymax),
             background.panel = "white",
             background.title = "grey95", 
             fontcolor.title  = "grey30",
             col.main         = "darkgray",
             col.title        = "darkblue",
             col.axis         = 1,
             col.baseline     = 1,
             title.width      = 1,
             cex.main         = 1,
             cex.title        = 0.5,
             cex.axis         = 0.4,
             baseline         = 0,
             innerMargin      = 0,
             frame            = TRUE)
}

### ----parallel-----------------------------------------------------------
if(image_fmt %in% c("png", "tiff")) {
  library(parallel)
  library(foreach)
  library(doParallel)
  cores <- 10
  clu   <- makeCluster(cores)
  registerDoParallel(clu, cores = cores)
  getDoParName()
  ptime <- system.time({    
    da  <- foreach (i = 1:nrow(sR_df), .packages = c("Gviz")) %dopar% { # For multiple CPUs
      options(ucscChromosomeNames=FALSE)
      pos       <- sR_df[i, ]
      #print(paste0(Sys.time(), " : ", pos$id))
      if(image_fmt == "png") {
        png_file <- paste0(outdir, '/', pos$id, '.png')
        png(png_file, width = 1600, height = 1200, res = 300)
      }else if(image_fmt == "tiff") {
        tiff_file <- paste0(outdir, '/', pos$id, '.tiff')
        tiff(tiff_file, width = 1600, height = 1200, res = 300, compression = "zip+p")
      }
      pts(gfrom = pos$start, gto = pos$end, main = pos$id)
      dev.off()
    }
    stopCluster(clu)
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}

###
if(image_fmt == "pdf"){
  ptime <- system.time({
    pdf_file <- paste0(outdir, '/', "seq_gviz_plots.pdf")
    pdf(pdf_file, width = 6, height = 5)
    for (i in 1:nrow(sR_df)) {
      options(ucscChromosomeNames = FALSE)
      pos <- sR_df[i, ]
      sd <- GRanges(seqnames = pos$seqnames, ranges = IRanges(start = pos$start, end = pos$end))
      ylimits <- mapply(function(x) {subcov <- subsetByOverlaps(cov[[x]], sd); max(values(subcov)$score)}, 1:8)
      ymax <- max(ylimits) * 1.2
      pts(gfrom = pos$start, gto = pos$end, ymax = ymax, main = pos$id)
    }
    dev.off()
  })
  print(paste0("Finish (create pngs): ", ptime[3]))
}
