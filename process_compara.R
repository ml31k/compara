################################################################
################################################################
## Name: process_compara.R
## Authors: Alejandra Medina <medina.alexiel@gmail.com>
##          "Mike" Minggao Liang <m.liang@mail.utoronto.ca>
## Date: January 22, 2014




################
##

## load libraries
library(GenomicRanges)
#library(gplots)
#library(gdata)
#library(gmodels)
#library(gtools)


################################################################
## processComparaData
##
## Function to process the list of input files that were generated using compara_match_peaks.pl
## compara_match.pl generates a table reporting per peak, the correspongin orthologous region in other organism
## according to the compara emsembl database annotation.
##
## Files with compara regions have in the first column the peak in the reference organims and then one column per queried organism
##
## One organims can have serveral regions matching to one peak in the reference organism, this can be caused by two reasons:
##  1) Compara alignment is split because is easier for the database to handel them this way, in this case the regions can be stitched back into one because they overlap.
##  2) There are paralogous regions in the genome, so you get two places where the sequence can align. After stitching the regions,
##     there are still several parts that are not cointiguos
##
## The function correct2 will stich together the aligment regions that overlap, in case of having several regions mapped, it return the longest one.
## The longest region has a better chance to overlap with a peak in an other organims.
##
## **Retrieving Paralogs**: In the future we would like to store and remmeber parelogs,
## so even if we will pick one for the comparison, we want to know which are the other regions

## inputFiles= Vector with the name of the files containing the output from compara_match_peaks.pl
## outputFile= Prefix for the outfile, this will be used to creat an Rdata object with a list containing ALL compara corrected files.
## org = organims list with short names. DEPRECATED, we will remove it.
## dborganisms = List of names of organims used in the analysis

do.paralogs <<- FALSE;

processComparaData <- function (fileNames, outputFile, dborganisms) {
    
################
## Declare global data.frame to contain the paralog regions 

    paralogs <<- data.frame(used.peak=NA, paralogs=NA)
  ##### Read in list of Compara input files) #####
    print (paste("Reading Compara Input Files"))
    print (fileNames)
#    fileNames = read.table(inputFiles, header=FALSE, as.is=TRUE, col.names="files", sep="\n", strip.white=FALSE) 

    ## Check the number of organims is the same as compara file, it has to match
    if (!(length(dborganisms)==length(fileNames)) ){
        print ("Error Number of organims is not the same as number of provided compara files")
        stop()
    }
    
 
  ##### Correcting the compara tables  ####
  # The program assumes you have run through the compara_match_peaks.pl script for EACH species set of ChIP data.
  # Output from compara_match_peaks.pl includes a header describing content.
    
    compara.list <- list()
    rejected.list <- list()
    for (i in 1:length(dborganisms)) {
        
    ## Upon importing the data it will apply the "correct" function to the file to format it correctly
    ## End result only allows one alignment per species per block
        print (paste("Reading compara file for",i,": ",fileNames[i])) 
        speciesTable <- read.table(fileNames[i],sep="\t",header=TRUE, comment.char = "", fill = TRUE) ## Read one compara file

    ## Correction will be done per organims, applaying function per column in the file 
        if(!grepl(dborganisms[i],speciesTable[1,1])){
            print ("Compara file does not correspond to the correct organism")
            print ("Files have to be provided in a consistent order of organims")
            print (paste("File", fileNames[i] ))
            print (paste("Organism" ,dborganisms[i] ))
            print (paste("Correct Orther:",dborganisms))
            stop()
        }
        print (paste("Correcting orthologous regions for:", fileNames[i] ))
        ## Apply correcting fucntion to stitch regions
        ##print (dim(speciesTable[,-1]))
        #cnames <- colnames(speciesTable[,-1])
        corrected.list <- apply(speciesTable[,-1], c(1,2), correct2.1) ## Correct per organism (column)
        #print (class(corrected.list))
        #print (dim(corrected.list))
        #print (head(corrected.list))

        ## Save the corrected table in a list with the organims name as ID
        compara.list[[dborganisms[i]]] =  corrected.list
    }

    ## Once corrections has been applied to all organims files the list with the information is saved in a Rdata object.
    save(compara.list, file=paste(outputFile, "-compara.list.RData", sep=""))
  
  
    ### saving compara object into separate tables for each species
    for (i in 1:length(dborganisms)) {
        print (paste("Printing compara list table for ",dborganisms[i]))
        #print (head(compara.list[[dborganisms[i]]] ))
        tempList = as.data.frame(compara.list[[dborganisms[i]]])
        write.table(tempList, file=paste(outputFile, "_compara.list.", dborganisms[i], ".txt", sep=""), append=FALSE, sep="\t", eol="\n", col.names=TRUE, quote=FALSE)  
    }
    
  return (compara.list)
}

################
## correct2
##
## One organims can have serveral regions matching to one peak in the reference organism, this can be caused by two reasons:
##  1) Compara alignment is split because is easier for the database to handle them this way, in this case the regions can be stitched back into one because they overlap.
##  2) There are paralogous regions in the genome, so you get two places where the sequence can align. After stitching the regions,
##     there are still several parts that are not cointiguos
##
## The function correct2 will stich together the aligment regions that overlap, in case of having several regions mapped, it return the longest one.
## The longest region has a better chance to overlap with a peak in an other organims.
##
## **Retrieving Paralogs**: In the future we would like to store and remmeber parelogs,
## so even if we will pick one for the comparison, we want to know which are the other regions

correct2 <- function (x, dborganism, stitch_bp=2){
    column.one.org <- x ## All matched peaks from the reference organims in one other organism
    column.one.org.correct <- unlist(lapply(column.one.org, correct2.1,  stitch_bp= stitch_bp  )) ## Correct per matched region
    return(column.one.org.correct)
}

## correct2.1
## Auxiliar function for correct2.
## One "lifted" peak will be splitted if the original peak in the reference genome matched more than one region in an other organism
## in this case severan regions will be separated by ":", each region will be stored as an independent one and then the group will be
## compared to see if they overlap, if they overlap (stitch_bp) they will be merge into one region.
## If after stitching the region there are still more than one ranges, only the longest one will be reported.
## The other one(s) are probably paralogs that in the future we would like to retrieve.

correct2.1 <- function(one.orth.peak,  stitch_bp=2){
    ## If region is NA return NA, this means there is no ortholog region for one peak.
    if (is.na(one.orth.peak)){
        return (NA)
        
    }  ## If the peak matched to several region in the other organism those will be separated by ":"
      ## These are the regions that will be "corrected"
    else if (grepl(":",one.orth.peak)){
        
        ## Separate each match 
        all.matches <- unlist(strsplit(one.orth.peak, split=":"))  
        split.matches.l <- lapply(lapply(all.matches,strsplit,split="\\."), unlist)
        split.matches <- do.call(rbind, split.matches.l)
        
        ## Store independent matches in one Genome Ranges object
        ortho.GR <- GRanges(seqnames = split.matches[,2], ranges = IRanges(as.numeric(split.matches[,3]),as.numeric(split.matches[,4])) )
        
        ## Stitch overlaping regions selected from the multiple sequence aligment that overlap by stitch_bp
        ortho.GR <- reduce(ortho.GR,min.gapwidth= stitch_bp)
        
        ## If two regions match the peak, select the longest one to be return
        ## In the future we want to remmember the other regions and report them
        if(length(ortho.GR)>1){
            longest.ortho.GR <- ortho.GR[which (width(ortho.GR)%in% max(width(ortho.GR)))]
            longest.ortho <- as.data.frame(longest.ortho.GR[1,])
            longest.ortho <- paste(split.matches[1,1],longest.ortho$seqnames,longest.ortho$start,longest.ortho$end, sep=".")     
            oneset.paralogs <- as.data.frame(ortho.GR)
            oneset.paralogs2 <-paste(split.matches[1,1],oneset.paralogs$seqnames,oneset.paralogs$start,oneset.paralogs$end, collapse=":",sep=".")
                                       # print (longest.ortho)
                                       # print (oneset.paralogs)
                                       # print (oneset.paralogs2 )
                                       # stop("DIED HERE")
            paralogs <<-rbind(paralogs,c(longest.ortho, oneset.paralogs2))
#            print (head(paralogs))
#            stop("DIED HERE")
        }
        
        else{
            longest.ortho <- as.data.frame(ortho.GR)
            longest.ortho <- paste(split.matches[1,1],longest.ortho$seqnames,longest.ortho$start,longest.ortho$end, sep=".")
            
        }
        
        return(longest.ortho) 
    }
    else{
        
        ## retunr region that didn't require correction
        return (one.orth.peak)
    }


}


 

################
## Function to find the peaks conserved between organisms
## peakFile - GRangesList object stored in a Rdata fila that will be loaded.
##            The list in fonformed by the independent Granges objects of the
##            peaks found with ChIP-seq in each organism
## compara.list - List containing the regions retrived from the compara alignment this region were "corrected" with the function correct2.
## org - organims ## DEPRECATED, consider removing
## dborganims - vector with the names of the organisms used throughout this analysis
## maxGap and MonOver - Passed to findOverlaps:  Intervals with a separation of ‘maxgap’ or less and
##          a minimum of ‘minoverlap’ overlapping positions, allowing for
##         ‘maxgap’, are considered to be overlapping. (Note from ?findOverlaps)
## outputFile - Prefix for the files to be created.
##            Tables containing the conservation information for each peak (Boolean variables),
##            and the corresponding peak in the other organims (one column per organism) that were considered the conserved region.
##            The column for the ref.organims will be filled with ones, DEPRECATED, remove after tool is stable 

analyseComparaData <- function(peakFile, compara.list, dborganisms, maxGap=0, minOver=1, outputFile) {
    load(peakFile) # Load GRanges List with original peaks
    print(paste("Ananlysing Compara data for organims:",dborganisms))
   
    ## Anlyse peaks for each organim (ref.org)
    for (i in 1:length(dborganisms)){
        ref.org <- dborganisms[i]
        print (paste("Comparing", ref.org ,"peaks with:"))
        local.ref.org.peaks.range <- peaks.GR[[ref.org]]
        ## Loop through all organisms (other.org) to compare the ref.org
        for (other.org in dborganisms){
            print (paste("    ","matched peaks from",other.org))
            if (ref.org == other.org){
                #By definition peaks in the reference organims are conserved in the reference organism
                mcols(local.ref.org.peaks.range)[other.org] <- rep("1", length(local.ref.org.peaks.range))
            }

            ## Declare vector that will contain the conservation information
            ## Vector is initalized in 0, so by defult we consider the peak will not be conserved in the other.org
            other.org.peaks.to.save <- rep("0", length(local.ref.org.peaks.range))
            other.org.peaks.to.save[is.na(compara.list[[ref.org]][,other.org])] = "X";
            other.org.peaks.paralog <- rep(NA, length(local.ref.org.peaks.range))

            ## Retrive for other.org the corresponding corrected compara table from the compara.list object
            other.org.in.ref.org <- compara.list[[other.org]][,ref.org]
            ## Get the corresponding column for the ref.org to be compared.
            ## Separeate the information into independent columns to be used in a GR object
            other.org.in.ref.org.l  <- lapply(lapply(other.org.in.ref.org,strsplit,split="\\."), unlist)
            ## Recover the original peak from the other.org
            other.org.in.ref.org.split <-do.call(rbind, other.org.in.ref.org.l)
            ##paste informaton together in one object
            other.org.in.ref.org.split <-cbind( other.org.in.ref.org.split , compara.list[[other.org]][,other.org])
            ## Remove NA comming from compara dabase when the region in other.org is not located in the multiplesequence aligment
            other.org.in.ref.org.split <- na.omit(other.org.in.ref.org.split)

            ## Creat genome Ranges object with the peak coordinates for the other.org "lifted to" ref.org 
            other.org.peaks.in.ref.GR <- GRanges(seqnames = other.org.in.ref.org.split[,2], ranges = IRanges(as.numeric(other.org.in.ref.org.split[,3]),as.numeric(other.org.in.ref.org.split[,4])), original.other.org.peak= other.org.in.ref.org.split[,5])

            ## Overlap ref.org peak ranges to the genome range object containing the other.org peaks in the ref.org 
            overlaps <- as.matrix(findOverlaps(local.ref.org.peaks.range, other.org.peaks.in.ref.GR, select="all", maxgap=maxGap, minoverlap=minOver))
           
            ## Loop through the peaks from the ref.org that overlapped to the "lift-over" peaks in the other.org,
            ## and retrieve the original peak in the other.org
            for (peak in unique(overlaps[,1])){
                other.peak <- overlaps[which( overlaps[,1] == peak),2]
                compara <- paste(mcols(other.org.peaks.in.ref.GR[other.peak])$original.other.org.peak,  sep=":") # Store corresponding peak in the other organism
                ## Recover paralogs from the paralog file table
#                print (compara)

                if (do.paralogs && grepl(compara, paralogs)){
                    other.org.peaks.paralog[peak] <- paralogs$paralogs[grep(compara, paralogs$used.peak)]
                }
                other.org.peaks.to.save[peak] <- compara
            }
            
            
            ## Save conservation information of each peak in the ref.org in the peak GR object under the other.org metadata          
            mcols(local.ref.org.peaks.range)[other.org] <-  other.org.peaks.to.save
            original.r <- as.data.frame(mcols(local.ref.org.peaks.range)$rejected.ranges )
#            print (other.org.peaks.paralog)
#            print (ref.org)
#            print (other.org)
            other.org.peaks.paralog <- apply(cbind(original.r, other.org.peaks.paralog),1,paste,collapse=";")
            mcols(local.ref.org.peaks.range)["rejected.ranges"] <-other.org.peaks.paralog
                            
        }
       
        all.orgs.cons <- as.data.frame(mcols(local.ref.org.peaks.range)[,dborganisms])
        all.orgs.cons.flat <- apply( all.orgs.cons,2, flat.conservation) ## flat, transform conservation information into a boolean value 1, conserved, 0 not conserved
        
        cons <- apply(all.orgs.cons.flat,1,paste,collapse="", sep="") ## Collapse conservation information into a boolean variable with as many posstions as organims,
                                                                      ## sorted as the organims are sorted
        
        mcols(local.ref.org.peaks.range)["conservation"] <- cons  ## assing conservation vector to metadata to be printed
        out.table.file <- paste(outputFile, "_peakConservation_", ref.org, ".txt", sep="")  ## Name of the outfile with the conservation information for the ref.org
        print(paste("Conservation Results for:" ,ref.org, "  ", out.table.file))
       
        write.table(as.data.frame(local.ref.org.peaks.range), file=out.table.file, append=FALSE, sep="\t", eol="\n", col.names=TRUE, row.names=FALSE, quote=FALSE)
        peaks.GR[[ref.org]] <- local.ref.org.peaks.range
        
        rm( local.ref.org.peaks.range)
    }
}

## Function called by analyseComparaData
## Chance all conservation annotation into a boolean variable
## Conservation originaly contains the name of the peak in the other.organims
## the function will change this into 1 for conserved keeping the 0 for the non-conserved
flat.conservation <- function(x){
    x[which(x!=0 & x!="X")] <- 1
    return(x)
}

