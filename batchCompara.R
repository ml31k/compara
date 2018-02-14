################################################################
################################################################
## Name: process_compara.R
## Authors: Kalvin Monk
##          Alejandra Medina <medina.alexiel@gmail.com>
##          "Mike" Minggao Liang <m.liang@mail.utoronto.ca>
## Date: January 22, 2014





## Batch compiler for running compara data experiments
## This program will run all the required steps to analyse the conservation of a set of ChIP-seq peaks
## The analysis can be re-run numerous times with different overlap/gap variables


## processComparaData which builds a compara.list file. This reformatting is the longest part of the script
## analyseComparaData which looks at the overlaps and generates the conservation data

################
## Load required functions

## Functions to process peak files in GFF format, output is a GenomeRangesList containing all GenomeRanges
#source ("peakConverter.R")

## Functions to process the outpur form compara_match_peaks.pl 
#source ("process_compara.R")


################
## batchComparaData
## inputInstructions - Tab delimited file 1)file containing the peak files list 2)file containing the compara_out file list
## 3)output prefix 4)maxOverlap:MinOverlap to match conserved peaks.


batchComparaData <- function (inputInstructions)
{
    ## Read table with input instructions
    inputFiles = read.table(inputInstructions, header=FALSE, as.is=TRUE, sep="", strip.white=FALSE)
    print(inputFiles)
    names(inputFiles) = c("fileList","outputFormat","variables")
    
    finalPeaks = list(NA)
    print("loaded")

    ## Several sets of runs can be contained in the inputInstructions file
    ## each with a list of peak files and  list of compara files
    for (i in 1:nrow(inputFiles)) {
        dataFiles = read.table(inputFiles$fileList[i], header=FALSE, as.is=TRUE, sep="", strip.white=FALSE)
        names(dataFiles) = c("orgs", "comparaFiles", "peakFiles")
        varList = processVariables(inputFiles$variables[i])
        dborganisms = dataFiles$orgs
        
        ## Process one set of peaks stored in one input file  
        peaks.GR = gffToPeakfile(dataFiles$peakFiles,inputFiles$outputFormat[i],organisms=dborganisms)
        print(paste("GR file generated with files in:", dataFiles$peakFiles))
        
        ## Process one set of outfiles form compara_match_peaks.pl in one input file 
        compara.list = processComparaData(dataFiles$comparaFiles, 
            inputFiles$outputFormat[i], 
            dborganisms)
        print(paste("Compara data processed with files in:",dataFiles$comparaFiles))
        
        ## For each combination of maxOVerlaps and MinOVerlaps values
        ## Run anayseComparaData
        ## Keep track of which combos are already done so as to avoid uneccessary runs

        done <-list("NA");
      
        for (j in 1:length(varList)) {
          
            maxGap <- as.integer(varList[[j]][1])
            minOver <- as.integer(varList[[j]][2])
            output <- inputFiles$outputFormat[i]
            mmID <- paste(output, "max", maxGap, "min", minOver, sep="")
            ## Checks if this combination has already been done, if so, skip.
            if (grepl(mmID, done)) {
                next()
            }
            done <- paste(done, mmID, sep="")
            ## Run Analyse compara with the corresponding data set and values 
            analyseComparaData(paste(output, "-peakFile.RData", sep=""),
                                        # peaks.GR, 
                               compara.list, 
                               dborganisms, 
                               maxGap, minOver, 
                               mmID                             
                               )
        }
    }
    print("finished!")

#    print(paralogs)
}

##### processVariables #####
# A quick function to split the maxgap and minoverlap variables into pairs for later use

processVariables = function(vars) {
  
  gapOver = unlist(strsplit(vars, ";"))
  gapOver = strsplit(gapOver, ":")
  return(gapOver)
  
}
