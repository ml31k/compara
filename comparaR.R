library(GenomicRanges)

coord2GR <- function(x){
    split.coord <- unlist(strsplit(x,"\\."))
    GRanges(seqnames=split.coord[2], ranges=IRanges(start=as.numeric(split.coord[3]), as.numeric(split.coord[4])))
}

stitchbp = 2 ## Compara often reports regions with small gaps. Regions separated by <= stitchbp will be merged

args <- commandArgs(trailingOnly=TRUE)

#args <- "~/mdwilson/mliang/projects/ChIP_exo_mm10_20160430/data/GEM_peaks/GEMbed/compara/CEBPA/rscript_instructions.txt"

if(length(args)==0){
    stop("No instruction file supplied. File should be in the form\npath_to_inputs\tTF_name\tmaxSep:minOverlap")
} else {
    instruction.file <- as.character(read.table(args[1],sep="\t",stringsAsFactors=FALSE)[1,])
    matchType <- args[2]
}

print(args)

if(is.null(matchType) | is.na(matchType)){
    matchType <- "P"
    print("No matchType specified: Matching in (P)ermissive mode")
} else if (matchType == "L"){
    print("matchType specified as 'L': Matching in (L)egacy mode")
} else if (matchType == "P"){
    print("matchType specified as 'P': Matching in (P)ermissive mode")
} else if (matchType == "S"){
    print("matchType specified as 'S': Matching in (S)trict mode")
} else {
    stop("unrecognized matchType")
}


file.table <- read.table(instruction.file[1], row.names=1, stringsAsFactors=FALSE)
colnames(file.table) <- c("compara_out","original_peak")

tf <- instruction.file[2]

maxgap <- as.numeric(gsub(":.*","",instruction.file[3]))
minovl <- as.numeric(gsub(".*:","",instruction.file[3]))

species.list <- row.names(file.table)

list.of.compara <- lapply(species.list, function(a.species){
    print(paste("Processing Compara Regions for:", a.species))
    a.compara.tab <- read.table(file.table[a.species,]$compara_out, comment.char="",header=T, stringsAsFactors=FALSE)
    a.compara.GR.tab <- lapply(species.list, function(other.species){
        print(other.species)
        other.species.regions <- a.compara.tab[,other.species]
        names(other.species.regions) <- a.compara.tab$X.Original_peak
        lifted.peaks <- split(rep(GRanges(seqnames="chr0",ranges=IRanges(1,1)), #oripk=NA),
                                  nrow(a.compara.tab)),1:nrow(a.compara.tab))
        single.regions <- t(sapply(other.species.regions[which(!grepl(":", other.species.regions) &
                                                             !is.na(other.species.regions))], function(x){
                                                                 unlist(strsplit(x, "\\."))[2:4]
                                                             }))
        single.ranges <- GRanges(seqnames=single.regions[,1],
                                 ranges=IRanges(as.numeric(single.regions[,2]),as.numeric(single.regions[,3])))                                                             
        lifted.peaks[which(!grepl(":", other.species.regions) & !is.na(other.species.regions))] <- 
            split(single.ranges, 1:length(single.ranges))
        multi.regions <- do.call(rbind,
                                 lapply(names(other.species.regions)[which(grepl(":", other.species.regions))],
                                        function(a.name){
                                            a.region <- other.species.regions[a.name]
                                            out <- matrix(unlist(strsplit(a.region,"[\\.:]")),ncol=4, byrow=T)
                                            out <- cbind(out, rep(a.name, nrow(out)))
                                        }))
        multi.ranges <- GRanges(seqnames=multi.regions[,2], oripk = multi.regions[,5],
                                ranges=IRanges(as.numeric(multi.regions[,3]), as.numeric(multi.regions[,4])))
        lifted.peaks[which(grepl(":", other.species.regions))] <-
            multi.reduced <- reduce(split(multi.ranges, multi.ranges$oripk), min.gapwidth=stitchbp)

            
        names(lifted.peaks) <- a.compara.tab$X.Original_peak
        return(lifted.peaks)        
    })
    names(a.compara.GR.tab) <- species.list
    return(a.compara.GR.tab)
})

names(list.of.compara) <- species.list

original.peaks <- lapply(species.list, function(a.species){
    peaks <- read.table(file.table[a.species, "original_peak"])
    peak.ranges <- with(peaks, GRanges(seqnames=V1, ranges=IRanges(V2, V3),
                                       oripk = paste(V1,V2,V3, sep=".")))
})
names(original.peaks) <- species.list

matched.compara.list <- lapply(species.list, function(origin.species){
    matched.matrix <- sapply(species.list, function(other.species){
        print(paste("Matching peaks between",origin.species, other.species))
        ori.lift2.other.ol <- as.data.frame(findOverlaps(list.of.compara[[origin.species]][[other.species]],
                                                         original.peaks[[other.species]], maxgap=maxgap, minoverlap=minovl))
        ori.lift2.other.ol$oripk <- paste0(other.species, ".",
                                           original.peaks[[other.species]]$oripk[ori.lift2.other.ol$subjectHits])

        other.lift2.ori.ol <- as.data.frame(findOverlaps(original.peaks[[origin.species]],
                                                         list.of.compara[[other.species]][[origin.species]],
                                                         maxgap=maxgap, minoverlap=minovl))
        other.lift2.ori.ol$oripk <- names(list.of.compara[[other.species]][[origin.species]])[other.lift2.ori.ol$subjectHits]

        if(matchType == "S"){
            conserved.pks <- intersect(ori.lift2.other.ol$queryHits, other.lift2.ori.ol$queryHits)
        } else if (matchType == "P") {
            conserved.pks <- unique(c(ori.lift2.other.ol$queryHits, other.lift2.ori.ol$queryHits))
        } else if (matchType == "L") {
            conserved.pks <- other.lift2.ori.ol$queryHits
        }

        merged.ol <- unique(rbind(ori.lift2.other.ol, other.lift2.ori.ol))
        
        out <- rep("X", length(original.peaks[[origin.species]]))
        out[sum(width(list.of.compara[[origin.species]][[other.species]]) > 1) >=1] <- "0"
        out[conserved.pks] <- sapply(conserved.pks, function(a.pk){
            paste(unique(subset(merged.ol, queryHits==a.pk)$oripk), collapse=":")
        })
        return(out)
    })
})
names(matched.compara.list) <- species.list

compara.out.list <- lapply(species.list, function(a.species){
    out <- as.data.frame(original.peaks[[a.species]])[,1:5]
    out$conservation <- apply(matched.compara.list[[a.species]],1, function(a.pk) paste0(gsub(".*\\..*","1",a.pk), collapse=""))
    out$rejected.ranges <- NA
    return(cbind(out, matched.compara.list[[a.species]]))
})

for(one.species in species.list){
    write.table(x=compara.out.list[[one.species]],
                file=paste0(tf,"max",maxgap,"min",minovl,"_peakConservation",one.species,".txt"),
                row.names=F, quote=F, sep="\t")
}
