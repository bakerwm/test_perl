## ~40s for each seq
## solution for big data:
## Skip the first N rows, (skip=N), read N rows, (nrow=N)
##

#getCovData <- function(CovDir="Rv.coverage", Strain="H37Rv", ID, Begin=1, End=10, Strand="+") {
getCovData <- function(CovDir="Rv.coverage", Input.Info) {
    CovDir <- "Rv.coverage"
    #Input.Info <- reads[1, ]
    ID     <- as.character(Input.Info$ID)
    Strain <- as.character(Input.Info$Strain)
    Begin  <- as.numeric(Input.Info$Begin)
    End    <- as.numeric(Input.Info$End)
    Strand <- as.character(Input.Info$Strand)
        
## Check the strand
    if(Strand == "+") {
        Strand <- ".p$"
    } else if(Strand == "-") {
        Strand <- ".n$"
    } else {
        stop("Strand should be +/-:")
    }
    CovFiles <- dir(CovDir, Strain)
    CovFiles <- CovFiles[grepl(Strand, CovFiles)] # Select the strand
    #Change the name to <StrainName>.<LibType>.<Strand>, eg: H37Rv.45SE.n
    CovFiles.trim <- sub("_(\\d+\\wE)\\.coverage", "\\.\\1", CovFiles, perl=TRUE)
    NameSplit <- strsplit(CovFiles.trim, "\\.")
    ##Create CovFiles Index
    CovFiles.index <- data.frame(matrix(unlist(NameSplit), length(CovFiles), byrow=TRUE))
    names(CovFiles.index) <- c("Strain", "LibType", "Strand")
    CovFiles.index$Path <- paste(CovDir, CovFiles, sep="/")
    CovFiles.index$Filename <- CovFiles.trim
    Lib.index <- data.frame(LibType=c("45SE", "81SE", "140PE", "200PE"), 
                            LibName=c("1. 18-40 nt", "2. 40-80 nt", "3. 80-140 nt", "4. >140 nt"))
    
    StartRow <- Begin - 50; if(StartRow < 0) StartRow <- 1 
    ReadRows <- End - Begin + 100
    
    ReadAllLibData <- data.frame(Position=NA, Coverage=NA, LibName=NA, ID=NA)
    for(m in 1:length(CovFiles.index$Path)) {
        Input.file <- CovFiles.index$Path[m]
        Temp.Data <- read.table(Input.file, skip=StartRow, nrows=ReadRows, comment.char="",
                               colClasses=c("character", "numeric", "numeric"))
        names(Temp.Data) <- c("Chr", "Position", "Coverage")
        FileData <- subset(Temp.Data, Position >= Begin & Position <= End,
                           select=c(Position, Coverage))
        LibName <- Lib.index$LibName[Lib.index$LibType == CovFiles.index$LibType[m] ]
        FileData$LibName <- LibName
        FileData$ID      <- ID
        ReadAllLibData <- rbind(ReadAllLibData, FileData)
    }
    ReadAllLibData <- na.omit(ReadAllLibData)
}

#################################
#################################
library("ggplot2")
tm <- theme(axis.line  = element_line(colour="black"),
            axis.title = element_text(size=16, colour="black"),
            axis.text  = element_text(size=12, colour="black"),
            plot.title = element_text(size=20, colour="black"),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(fill="gray70", colour="black", size=0.5) )

#reads <- read.table("test.txt", comment.char="")
reads <- read.table("CYRNAseqView.txt", comment.char="")
names(reads) <- c("ID", "Strain", "Length", "Begin", "End", "Strand")

pdf("Rplot.pdf", width=8, height=6)
system.time(
    for(i in 1:nrow(reads)) {
        Line.Info <- reads[i, ]
        Begin <- Line.Info$Begin; End <- Line.Info$End; Length <- Line.Info$Length
        Linedata <- getCovData(CovDir="Rv.coverage", Input.Info=Line.Info)  ## Function
        LinePosition <- paste(Line.Info$Begin, Line.Info$End, Line.Info$Strand, Line.Info$Length, sep=":")
        figtitle <- paste(Line.Info$Strain, Line.Info$ID, LinePosition, sep=" ")
        p1 <- ggplot(Linedata, aes(x=Position, y=Coverage)) + ggtitle(figtitle) + tm +
            geom_area(fill="black") + scale_x_continuous(breaks=seq(Begin, End, Length/7))
        print(p1 + facet_grid(LibName ~ ID) )
        print(p1 + facet_grid(LibName ~ ID, scales="free_y") )
    }
)
dev.off()

