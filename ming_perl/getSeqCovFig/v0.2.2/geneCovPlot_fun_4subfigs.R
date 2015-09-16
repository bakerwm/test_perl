### ----plot the sRNA on reference genome----------------------------------
### ----four libraries
### ----parse coverage file, index-----------------------------------------
getcovidx <- function(covdir, chr = "H37Rv", readbig = FALSE) {
  readbigfile <- function(x, bin, desc, covdir) {
    library(bigmemory)
    read.big.matrix(x, sep = "\t", header = FALSE, type = "integer",
                    col.names = c("chr", "position", "count"),
                    backingpath = dirname(bin),
                    backingfile = basename(bin),
                    descriptor  = basename(desc))
  }
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

### ----Generate plots for input seq---------------------------------------
cov2fig <- function(subcov, pos, extend = 200, yfree = FALSE){
  library(ggplot2)
  library(scales)
  ### transform lib name
  libname <- c("01" = "18-40 nt", "02" = "40-80 nt",
               "03" = "80-140 nt", "04" = ">140 nt")
  subcov$libtype <- as.factor(libname[subcov$libtype])
  subcov$libtype <- factor(subcov$libtype, levels = c("18-40 nt", "40-80 nt", "80-140 nt", ">140 nt"))
  #subcov$strand  <- as.factor(subcov$strand)
  #levels(subcov$strand)  <- factor(c("+", "-"))
  ### axix ticks, labels
  xstart  <- ceiling(min(subcov$position) / 10) * 10
  xend    <- floor(max(subcov$position) / 10) * 10
  xstep   <- floor(floor((xend - xstart) / 7) / 10) * 10
  if(xstep == 0) xstep <- 5
  xbreaks <- seq(xstart, xend, by = xstep)
  ylimits <- c(0, max(subcov$count) * 1.2)
  ybreaks <- round(seq(1, max(subcov$count), length.out = 3))
  ### for libraries
  libfigs <- lapply(levels(subcov$libtype), function(x) {
    libcov <- subset(subcov, libtype == x)
    if(yfree) {
      ylimits <- c(0, max(libcov$count) * 1.2)  
      ybreaks <- round(seq(1, max(libcov$count), length.out = 3))
    }
    mytheme <- theme_classic() +
      theme(plot.title        = element_text(size = 10, colour = "black",
                                             face = "bold"),
            #plot.margin       = unit(c(0.2, 0.2, 0.2, 0.4), "cm"),
            axis.title        = element_blank(),
            panel.border      = element_rect(size = 0.6, colour = "black", 
                                             fill = NA),
            legend.background = element_blank(),
            legend.text       = element_text(size = 8, colour = "black"),
            legend.title      = element_text(size = 10, colour = "black",
                                             face = "bold"),
            legend.position   = c(0.9, 0.75))
    p <- ggplot(libcov, aes(x = position, y = count, fill = strand, colour = strand)) +
      geom_area(position = "identity", alpha = 0.2) +
      geom_line(size = 0.5) 
    p1 <- p + geom_vline(xintercept = c(pos$start, pos$end), colour = "blue",
                         size = .2, linetype = "longdash") +
      scale_color_manual(values = c("red", "blue")) +
      scale_x_continuous(breaks = xbreaks) +
      scale_y_continuous(limits = ylimits,
                         breaks = ybreaks, 
                         labels = scientific_format(digits = 0)) +
      ggtitle(x) +
      mytheme
    
    return(p1)
  })
  return(libfigs)
}
###
pos2fig   <- function(pos, covidx, extend = 200, yfree = FALSE) {
  subcov <- data.frame()
  for(i in 1:nrow(covidx)) {
    ### define the range of the x-axis:
    start    <- pos$start - extend
    if(start < 0) start <- 1
    end      <- pos$end + extend
    df1      <- attach.big.matrix(covidx$desc[i], covidx$bin[i]) # read desc, bin
    df2      <- df1[mwhich(df1, c("position"), list(c(start, end)),
                           list(c('ge', 'le')), 'AND')]
    df2      <- as.data.frame(df2)
    df2$id   <- pos$id
    df2$chr  <- pos$chr
    df2$libtype <- covidx$libtype[i]
    df2$strand  <- covidx$strand[i]
    subcov <- rbind.fill(subcov, df2)
    rm(start, end, df1, df2)
  }
  ### pos
  figlibs <- cov2fig(subcov, pos, extend = extend, yfree = yfree)
  return(figlibs)
}

readpos <- function(x) {
  df <- read.table(x, header = FALSE, sep = "\t", na.strings = "NA")
  df2 <- df[, 1:6]
  names(df2) <- c('id', 'chr', 'length', 'start', 'end', 'strand')
  df2 <- within(df2,  note <- paste0(chr, '_', id, ' (', length, ' nt) ', start, ':', end, ':', strand))
}
