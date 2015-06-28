## This function is designed to pickout coverage data from 4 libraries of a position range.

## Function 1: read coverage data
## Funcioon 2: get subset of data for one read
getCovData <- function(Strain="H37Rv", CovDir="Rv.coverage") {
    #CovDir <- "Rv.coverage"; Strain <- "test"
    CovFiles <- dir(CovDir, Strain)
    #Change the name to <StrainName>.<LibType>.<Strand>, eg: H37Rv.45SE.n
    CovFiles.trim <- sub("_(\\d+\\wE)\\.coverage", "\\.\\1", CovFiles, perl=TRUE)
    NameSplit <- strsplit(CovFiles.trim, "\\.")
    ##Create CovFiles Index
    CovFiles.index <- data.frame(matrix(unlist(NameSplit), length(CovFiles), byrow=TRUE))
    names(CovFiles.index) <- c("Strain", "LibType", "Strand")
    CovFiles.index$Path <- paste(CovDir, CovFiles, sep="/")
    CovFiles.index$Filename <- CovFiles.trim
    
    ## Read Coverage files
    CovData <- list()
    for(i in 1:length(CovFiles)) {
        Filename <- CovFiles.index[i, ]$Filename
        FilePath <- CovFiles.index[i, ]$Path
        CovData[[Filename]] <- read.table(FilePath, comment.char="", 
                                          colClasses=c("character", "numeric", "numeric"))
        names(CovData[[Filename]]) <- c("Chr", "Position", "Coverage")
    }
    ## Create Library Index
    Lib.index <- data.frame(LibType=c("45SE", "81SE", "140PE", "200PE"), 
                            LibName=c("1. 18-40 nt", "2. 40-80 nt", "3. 80-140 nt", "4. >140 nt"))
    CovData[["Lib.index"]] <- Lib.index
    CovData[["Fileinfo"]]  <- CovFiles.index
    CovData
}

getReadCovData <- function(CovData, Input.Info){
    # Input.Info is a data.frame, with ID/Strain/Length/Begin/End/Strand columns
    Strand <- Input.Info$Strand
    Begin  <- Input.Info$Begin
    End    <- Input.Info$End
    ## Check the strand
    if(Strand == "+") {
        Strand <- ".p$"
    } else if(Strand == "-") {
        Strand <- ".n$"
    } else {
        stop("Strand should be +/-:")
    }
    LibChoose <- names(CovData)
    LibChoose <- LibChoose[grepl(Strand, LibChoose)]
    Fileinfo  <- CovData[["Fileinfo"]]
    Lib.index <- CovData[["Lib.index"]]
    ## Create cov data.frame of the read
    ReadAllLibData <- data.frame(Position=NA, Coverage=NA, LibName=NA, ID=NA)
    for(i in 1:length(LibChoose)) {
        ReadLibData <- subset(CovData[[LibChoose[i]]], Position>=Begin & Position<=End, 
                              select=c(Position, Coverage))
        LibType <- Fileinfo$LibType[Fileinfo$Filename == LibChoose[i]]
        LibName <- Lib.index$LibName[Lib.index$LibType == LibType]
        ReadLibData[["LibName"]] <- LibName
        ReadLibData[["ID"]] <- Input.Info$ID
        ReadAllLibData <- rbind(ReadAllLibData, ReadLibData)
    }
    ReadAllLibData <- na.omit(ReadAllLibData) 
}
