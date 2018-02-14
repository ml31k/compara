################################################################
## Take as input the compara output and spli it in separated
## bed files for conserved and not conserved regions

## The progam read in the column names and it will identify the
## anlysed organisms and creat the corresponding shared file per
## organism

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
compara.output.file <-args[1]
compara.prefix.bed <-args[2]
ref.org <- args[3]
print(compara.output.file)
print(compara.prefix.bed)

#compara.output.file <-"./results/compara70/20141121/20141103_mergedBioRep/AoEC_none/Rad21/Rad21max0min1_peakConservation_homo_sapiens.txt"
#compara.prefix.bed <-"./results/compara70/20141121/20141103_mergedBioRep/AoEC_none/Rad21/Rad21max0min1_peakConservation_homo_sapiens"

#compara.output.file <-"../results/compara70/20141121/20141103_mergedBioRep/AoEC_none/H3K36me3/H3K36me3max0min1_peakConservation_homo_sapiens.txt"
#compara.prefix.bed <-"../results/compara70/20141121/20141103_mergedBioRep/AoEC_none/H3K36me3/H3K36me3max0min1_peakConservation_homo_sapiens"
#compara.output.file <- "../results/compara70/20141121/20141103_mergedBioRep/AoEC_none/cJun/cJunmax0min1_peakConservation_homo_sapiens.txt"
#compara.prefix.bed <-  "../results/compara70/20141121/20141103_mergedBioRep/AoEC_none/cJun/cJunmax0min1_peakConservation_homo_sapiens"
#, colClasses =column.class
compara.output <- read.table(file=compara.output.file, stringsAsFactor=FALSE, header=TRUE)


## The basic information is in the first 7 column
## if the number is higher several animal are being analyzed
all.columns <- colnames(compara.output)

## get organisms
organisms <- all.columns[-c(1:7)]
num.organism <- length(organisms)

central.organism <- ref.org
##Get table in bed format
compara.output$seqnames <- paste("chr",sep="",compara.output$seqnames)
compara.output$strand <- "+"

## bed file name for all peaks
## print all peaks in bed format

all.peaks.bed.file <- paste(compara.prefix.bed,"_allpeaks.bed", sep="")
write.table(compara.output[,c(1,2,3,6,4,5)],file=all.peaks.bed.file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## if analysis is only in one organism stop
if (num.organism==1){
    print ("Only one organisms available in the analysis no conservation was reported")
    stop()
}


## get peaks conserved in any organisms
other.org <- organisms[!(organisms==central.organism)]

print (paste ("other organisms:", paste(other.org)))
print (length(other.org))
if(length(other.org)==1){
    conserved.any <- compara.output[!(compara.output[,other.org]=="X" | compara.output[,other.org]=="0") ,]
    conserved.notcons <- compara.output[(compara.output[,other.org]=="X" | compara.output[,other.org]=="0"),]
}else {
    conserved.any <- compara.output[apply(compara.output[,other.org],1,function(x) any((x!= "X")& (x!="0"))),]
    conserved.notcons <- compara.output[!apply(compara.output[,other.org],1,function(x) any((x!= "X")& (x!="0"))),]

}
anyconserved.bed.file <- paste(compara.prefix.bed,"_ConservedAny.bed", sep="")
write.table(conserved.any[,c(1,2,3,6,4,5)],file=anyconserved.bed.file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## get not conserved peaks
notconserved.bed.file <- paste(compara.prefix.bed,"_NotConserved.bed", sep="")
write.table(conserved.notcons[,c(1,2,3,6,4,5)],file=notconserved.bed.file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## Check point

if ((dim(conserved.any)[1]+dim(conserved.notcons)[1])!=dim(compara.output)[1]){
    stop("checkpoint failed")
}

## loop organisms conservation

for (one.org in other.org){
    print(one.org)
    
    conserved.one.org <- compara.output[compara.output[,one.org]!= "X"&compara.output[,one.org]!= "0",]

    org.file <- paste(compara.prefix.bed,"_Conserved_",one.org,".bed", sep="")
    
#    if (one.org == "bos_taurus"){
#        
#        org.file <-  paste(compara.prefix.bed,"_ConservedCow.bed", sep="") 
#        
#    }else if(one.org == "rattus_norvegicus"){
#        org.file <-  paste(compara.prefix.bed,"_ConservedRat.bed", sep="") 
#        
#    }else if(one.org == "homo_sapiens"){
#        org.file <-  paste(compara.prefix.bed,"_ConservedHuman.bed", sep="") 
#    }
    
    print(org.file)
    write.table(conserved.one.org[,c(1,2,3,6,4,5)],file=org.file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}

## Get conserved in all
if (length(other.org)>1){
    cons.all <- paste(rep(1,length(other.org)), collapse="")
    compara.cons.all.org <- compara.output[grep(paste(cons.all,"$", sep=""), compara.output$conservation),]
    org.all.file <-  paste(compara.prefix.bed,"_ConservedAll.bed", sep="")
    print(org.all.file)
    write.table(compara.cons.all.org[,c(1,2,3,6,4,5)],file=org.all.file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}
