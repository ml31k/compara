################################################################
################################################################
## Name: process_compara.R
## Authors: Alejandra Medina <medina.alexiel@gmail.com>
##          "Mike" Minggao Liang <m.liang@mail.utoronto.ca>
## Date: January 22, 2014


################
## Load Genomic Ranges
library(GenomicRanges)

################
## gffToPeakfile
##
## A list of known ChIP-seq datapoints pulled in from the GFF files from each
## The goal is to import each species into a separate GenomicRange and then amalgamate those into a genomic range list
## Of the format: chromosome, source, type, start, stop, unknown, strand, unknown, name
## Organims should be provided in the same order as their corresponding files
##
## peakFileNames = file containing the names of the peak files for each organims in the corresponding order.One line per file
## ourFile = Prefix to store the GenomeRangesList object with the peaks from all organisms
## organims = vector with the name of the organims in the corresponing orther
gffToPeakfile = function(fileNames, outputFile,organisms) {

  # First read the peakFileNames from the list. This will allow a little more flexibility/ease in running the program
  # If you are dealing with multiple ChIP-seq data types, you can supply a different set of file names for 
  # each type of Transcription Factor
  
#  fileNames = read.table(peakFileNames, header=FALSE, as.is=TRUE, col.names="files", sep="\n", strip.white=FALSE)
  
  csGR.list = list();
  print (organisms)
  if (!(length(organisms)==length(fileNames)) ){
      print ("Error Number of organims is not the same as number of provided peak files")
      stop()
  }
  
  for (i in 1: length(fileNames)) {
  
      ## Read each peak file, assigment to genome will be based on the organisms object
      print (paste("Reading peak file", fileNames[i],"for organism", organisms[i]) )

      ## Read file into a table
      csData <- read.table(fileNames[i], sep="", header=F)
      names(csData) <- c("chr", "mStart", "mEnd", "seqName")

      ## Transform coordinates into IRange object
      csRanges <- IRanges(start=csData$mStart, end=csData$mEnd, width=NULL, names=csData$seqName)
      cons <- data.frame(conservation=rep(NA,length(csRanges))) ## Declare empty conservation vector

      ## Declare genome ranges object for the peaks.
      csGR  <- GRanges(seqnames=csData$chr, ranges=csRanges, conservation=rep(NA, length(csRanges)), rejected.ranges=rep(NA, length(csRanges)))

      ## Delacre empty columns to store peak conservation in each organims
      empty.org.report <- as.data.frame(matrix(nrow = length(csRanges), ncol = length(organisms)))
      colnames(empty.org.report) <- organisms

      ## Attach the organims column to the Peak.GR 
      values(csGR) <-cbind(values(csGR),DataFrame(empty.org.report))

      ## Store the peak.GR in aList
      csGR.list [[organisms[i]]] <- csGR
  }
 ## Store the peak.GR in a GenomeRangesList
  peaks.GR <- GRangesList(csGR.list)

  ## Sabe the GenomicRangesList in a Rdata file
  file.name <- paste(outputFile, "-peakFile.RData", sep="")
  print (paste("Saving GenomeRanges file for peaks", file.name))
  save(peaks.GR, file=file.name)
  
  return(peaks.GR)
  
}
