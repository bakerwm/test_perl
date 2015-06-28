# Extract subset data, Using bigmemory package
getCovData <- function(FileIndex, InPos) {
    ## Check the strand
    subIndex <- subset(FileIndex, Strand == InPos$StrType)
    
    # Expand the range of input seq
    StartRow <- InPos$Begin - 200; if(StartRow < 0) StartRow <- 0
    EndRow   <- InPos$End + 200
    ReadRows <- InPos$Length + 400
    
    AllLibData <- data.frame(Chr      = NULL,
                             Position = NULL, 
                             Count    = NULL, 
                             LibName  = NULL, 
                             ID       = NULL)
    
    for(n in 1:nrow(subIndex)) {
        LibBigdesc <- dget(subIndex$Bigdesc[n])
        LibData <- attach.big.matrix(LibBigdesc)
        LibData <- LibData[mwhich(LibData, c("Position"), list(c(StartRow, EndRow)),
                                  list(c('ge', 'le')), 'AND')]
        LibData <- as.data.frame.matrix(LibData)
        LibData$LibName <- subIndex$LibName[n]
        LibData$ID      <- InPos$ID
        LibData$Chr     <- subIndex$Strain[1]
        AllLibData      <- rbind(AllLibData, LibData)
    }
    AllLibData
}

#################################
library("ggplot2")
library("plyr")
library("bigmemory")

# Load parameters
Usage <- "Usage: Rscript getSeqCovFig.R <InDir> <InFile> <OutPdf>"
args     <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop(Usage)
    
InDir    <- args[1]
InFile   <- args[2]
OutPdf   <- args[3]
print( paste(c("Input para:", args), sep=" ") )

################
FlankGap <- 200
#InDir  <- "Rv.coverage"
#InFile <- "test.txt"
################

# Read in seq info
InLists <- read.table(InFile, comment.char="")
names(InLists) <- c("ID", "Strain", "Length", "Begin", "End", "Strand")
InLists$StrType <- mapvalues(InLists$Strand, c("+", "-"), c("p", "n")) 
InLists$Note   <- apply(InLists[, c(1,3:6)], 1, paste, collapse=":")

# Create Cov Index
CovFiles <- dir(InDir, ".coverage.")
if(length(CovFiles) < 2 ) stop ("Check cov files in : InDir")
CovNames <- strsplit(CovFiles, "\\.")
CovIndex <- data.frame(matrix(unlist(CovNames), length(CovFiles), byrow=TRUE))
names(CovIndex)   <- c("Strain", "LibType", "Cov", "Strand")
CovIndex$LibName  <- mapvalues(CovIndex$LibType, 
                               c("01","02","03","04"), 
                               c("1. 18-40 nt", "2. 40-80 nt", "3. 80-140 nt", "4. >140 nt"))
CovIndex$FilePath <- paste(InDir, CovFiles, sep="/")
CovIndex$Filename <- CovFiles
CovIndex$Bigbin   <- paste(CovIndex$Filename, ".bin", sep="") 
CovIndex$Bigdesc  <- paste(CovIndex$Filename, ".desc", sep="")

CovIndex <- CovIndex[c("Strain", "Strand", "LibType", "LibName", "Filename", 
                       "FilePath", "Bigbin", "Bigdesc")]

# Create bin & desc files for each coverage file
#dir.create(paste(InDir, "desc", sep="/"), showWarnings=FALSE)
CheckDesc <- dir(path="./", pattern="*.desc")
if(length(CheckDesc) != nrow(CovIndex)) {
    for(i in 1:nrow(CovIndex)) {
        tFile <- CovIndex$FilePath[i]
        tbin  <- CovIndex$Bigbin[i]
        tdesc <- CovIndex$Bigdesc[i]
        tmp   <- read.big.matrix(tFile, sep="\t", header=FALSE, 
                                 type="integer",
                                 col.names=c("Chr", "Position", "Count"),
                                 backingfile=tbin,
                                 descriptor=tdesc)
    }    
}

#pdf("ReadsCov.pdf", width=8, height=6)

pdf(OutPdf, width=8, height=6)
tm <- theme(axis.line  = element_line(colour="black"),
            axis.title = element_text(size=16, colour="black"),
            axis.text  = element_text(size=12, colour="black"),
            plot.title = element_text(size=20, colour="black"),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(fill="gray70", colour="black", size=0.5) )

system.time(    
    for(i in 1:nrow(InLists)) {
        LinePos  <- InLists[i, ]
        FigBegin <- LinePos$Begin; FigEnd <- LinePos$End
        LinePos$Begin  <- LinePos$Begin - FlankGap; if(LinePos$Begin < 0) LinePos$Begin <- 1
        LinePos$End    <- LinePos$End + FlankGap
        LinePos$Length <- LinePos$End - LinePos$Begin + 1
        LineData  <- getCovData(FileIndex=CovIndex, InPos=LinePos)
        LineTitle <- paste(LinePos$Strain, LinePos$Note, sep=" ")
        # Select subset of data
        LineData <- subset(LineData, 
                           Position >= LinePos$Begin & Position <= LinePos$End)
        # For labels on x-axis
        xStart <- 10 * floor(LinePos$Begin/10)
        xEnd   <- 10 * floor(LinePos$End/10)
        xStep  <- 10 * floor(floor(LinePos$Length/7)/10)
        xBreak <- seq(xStart, xEnd, by=xStep)        
        p1 <- ggplot(LineData, aes(x=Position, y=Count)) +
            geom_area(fill="black") + 
            geom_vline(xintercept=c(FigBegin, FigEnd), colour="blue",
                       linetype="longdash", size=.3) +
            scale_x_continuous(breaks=xBreak) +
            ggtitle(LineTitle) + tm
        
        print(p1 + facet_grid(LibName ~ ID) )
        print(p1 + facet_grid(LibName ~ ID, scales="free_y") )
    }
)

dev.off()
