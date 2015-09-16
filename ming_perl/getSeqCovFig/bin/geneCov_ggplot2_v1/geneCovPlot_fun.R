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
### ----Generate plots for input seq---------------------------------------
pos2cov <- function(pos, covidx, strandness = TRUE, extend = 200) {
  idx2cov <- function(x, str = "+", extend = extend){
    subidx <- covidx[covidx$strand == str, ]
    subcov <- data.frame()
    for(i in 1:nrow(subidx)) {
      ### define the range of the x-axis:
      start    <- x$start - extend
      if(start < 0) start <- 1
      end      <- x$end + extend
      df1      <- attach.big.matrix(subidx$desc[i], subidx$bin[i]) # read desc, bin
      df2      <- df1[mwhich(df1, c("position"), list(c(start, end)),
                             list(c('ge', 'le')), 'AND')]
      df2      <- as.data.frame(df2)
      df2$id   <- x$id
      df2$chr  <- x$chr
      df2$libtype <- subidx$libtype[i]
      df2$strand  <- str
      #rm(start, end, df1)
      subcov <- rbind.fill(subcov, df2)
    }
    return(subcov)
  }
  strands <- c("+", "-")
  if(strandness) strands <- x$strand
  outcov <- lapply(strands, function(x) idx2cov(pos, str = x, extend = extend))
  return(outcov)
}
###
cov2fig <- function(subcov, extend = 200, fill_color = "black"){
  library(ggplot2)
  library(scales)
  ### axix ticks, labels
  xstart  <- ceiling(min(subcov$position) / 10) * 10
  xend    <- floor(max(subcov$position) / 10) * 10
  xstep   <- floor(floor((xend - xstart) / 7) / 10) * 10
  if(xstep == 0) xstep <- 5
  xbreaks <- seq(xstart, xend, by = xstep)
  ybreaks <- round(seq(1, max(subcov$count), length.out = 3))
  ### transform lib name
  libname <- c("01" = "18-40 nt", "02" = "40-80 nt",
               "03" = "80-140 nt", "04" = ">140 nt")
  subcov$libtype <- as.factor(libname[subcov$libtype])
  levels(subcov$libtype) <- factor(c(">140 nt", "80-140 nt", "40-80 nt", "18-40 nt"))
  ###
  p <- ggplot(subcov, aes(x = position, y = count)) +
    geom_area(fill = fill_color) +
    geom_vline(xintercept = c(pos$start, pos$end), colour = "blue",
               size = .2, linetype = "longdash") +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks, 
                       labels = scientific_format(digits = 0)) +
    ggtitle(pos$note)
  p1 <- p + facet_grid(libtype ~ .) 
  p2 <- p + facet_grid(libtype ~ ., scales = "free_y")
  subfigs <- list(p1, p2)
  return(subfigs)
}
###
pos2fig   <- function(pos, covidx, extend = 200) {
  ### on both strand
  covlist <- pos2cov(pos, covidx, extend  = extend, strandness = FALSE)
  ### on plus 
  figsp   <- cov2fig(covlist[[1]], extend = extend, fill_color = "orangered")
  figsn   <- cov2fig(covlist[[2]], extend = extend, fill_color = "royalblue1")
  return(c(figsp, figsn))
}

readpos <- function(x) {
  df <- read.table(x, header = FALSE, sep = "\t", na.strings = "NA")
  df2 <- df[, 1:6]
  names(df2) <- c('id', 'chr', 'length', 'start', 'end', 'strand')
  df2 <- within(df2,  note <- paste0(chr, '_', id, ' (', length, ' nt) ', start, ':', end, ':', strand))
}
