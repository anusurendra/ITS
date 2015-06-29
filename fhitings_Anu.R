#!/usr/bin/Rscript

#check to see if pacages are installed and install if needed
ch1<-require(reshape)
ch2<-require(plyr)
ch3<-require(gtools)

if(ch1){
  library(reshape)
}else{
  install.packages("reshape")
}

if(ch2){
  library(plyr)
}else{
  install.packages("plyr")
}


if(ch3){
  library(gtools)
}else{
  install.packages("gtools")
}

#set the significant digits
options(scipen=2000)

#read in arguments
args <- commandArgs(trailingOnly = TRUE)

#Intialize functions
#Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")


printHelpStatement<- function(error.string){
  statement <- "USAGE\nfhitings [-h] [-filedir filepath] [-taxafile filename] [-evalue_cutoff float_value] [-perc_identity_cutoff float_value] [-blastfiles filenames]\n
  DESCRIPTION\n-filedir <String>\npath to Blast results directory\nDefault = '.'\n
  -taxafile <String>\nPath to the taxanomy assignment CSV file\nDefault = './Hogan_db_mod.csv'\n
  -evalue_cutoff <Float>\nEvalue cutoff for blast results\nDefault = '0.001'\n 
  -perc_identity_cutoff <Float>\npercentage cutoff for taxanomic assignment resolution\nDefault = '0.08'\n
  -ucfiles <String>\nlist of mapped reads to OTUs files to process separated by ',' \nDefault = ''\n
  -blastfiles <String>\nlist of blast output files to process separated by ',' \nDefault = ''\n"
  
  statement <- ifelse(is.null(error.string),statement,paste(error.string,statement,sep="\n"))
  return(statement)
}


#Get variables
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))


if(length(which(argsDF[,1]=="h"))==1){
  writeLines(printHelpStatement(NULL))
  quit(save = "no", status = 1, runLast = FALSE)
}


#Initialize variables
FILEDIR <- ifelse(length(which(argsDF[,1]=="filedir"))==1,as.vector(argsDF[which(argsDF[,1]=="filedir"),2]),getwd())
TAXAFILE <- ifelse(length(which(argsDF[,1]=="taxafile"))==1,as.vector(argsDF[which(argsDF[,1]=="taxafile"),2]),paste(getwd(),"Hogan_db_mod.csv",sep="/"))
EVALUE_CUTOFF <- ifelse(length(which(argsDF[,1]=="evalue_cutoff"))==1,as.double(as.vector(argsDF[which(argsDF[,1]=="evalue_cutoff"),2])),0.001)
HIT_ABUNDANCE_CUTOFF <- ifelse(length(which(argsDF[,1]=="hit_abundance_cutoff"))==1,as.double(as.vector(argsDF[which(argsDF[,1]=="hit_abundance_cutoff"),2])),0.8)


if(grepl(",",argsDF[which(argsDF[,1]=="blastfiles"),2])==TRUE){
  BLAST_FILES <- ifelse(length(which(argsDF[,1]=="blastfiles"))==1,strsplit(as.vector(argsDF[which(argsDF[,1]=="blastfiles"),2]),","),NULL)
}else{
  BLAST_FILES <- ifelse(length(which(argsDF[,1]=="blastfiles"))==1,as.vector(argsDF[which(argsDF[,1]=="blastfiles"),2]),NULL)
}

print(BLAST_FILES)

if(grepl(",",argsDF[which(argsDF[,1]=="ucfiles"),2])==TRUE){
  UC_FILES <- ifelse(length(which(argsDF[,1]=="ucfiles"))==1,strsplit(as.vector(argsDF[which(argsDF[,1]=="ucfiles"),2]),","),NULL)
}else{
  UC_FILES <- ifelse(length(which(argsDF[,1]=="ucfiles"))==1,as.vector(argsDF[which(argsDF[,1]=="ucfiles"),2]),NULL)
}

print(UC_FILES)


taxanomic.csv.m<-as.data.frame(read.delim(file=TAXAFILE,
                                          header=T,as.is=T,row.names=NULL,fill=F,sep=",",
                                          quote="",comment.char="",blank.lines.skip=F,strip.white=T),stringsAsFactors=F)


bf.vec <- unlist(BLAST_FILES);
ucf.vec <- unlist(UC_FILES);

#print(BLAST_FILES)
for(i in c(1:length(bf.vec))){
  
  
  bf <- bf.vec[i]
  ucf <- ucf.vec[i]
  
  print(file.exists(paste(FILEDIR, bf, sep="/")))
  print(paste("Processing ",bf,sep=""))
  
  if(!file.exists(paste(FILEDIR, bf, sep="/")) || !file.exists(paste(FILEDIR, ucf, sep="/"))){
    writeLines(printHelpStatement("Ensure that the blast/uc file exists."))
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  print(paste("Reading Blast Output",sep=""))
  blast.results.m<-as.data.frame(read.delim(file=paste(FILEDIR, bf, sep="/"),
                                            header=F,as.is=T,row.names=NULL,fill=F,sep="\t",
                                            quote="",comment.char="",blank.lines.skip=F),stringsAsFactors=F)
  
  
  species.count.m <- NULL
  genus.count.m <- NULL
  family.count.m <- NULL
  order.count.m <- NULL
  subclass.count.m <- NULL
  class.count.m <- NULL
  subphylum.count.m <- NULL
  phylum.count.m <- NULL
  kingdom.count.m <- NULL
  
  taxanomic.assignment.m <- NULL
  
  
  colnames(blast.results.m) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  
  print(paste("Applying evalue filter",sep=""))
  pass.evalue.filter <- which(blast.results.m[,"evalue"]<=EVALUE_CUTOFF)
  
  blast.results.m.evalue.filter <- blast.results.m[pass.evalue.filter,]
  blast.results.m.evalue.filter.list <- split(blast.results.m.evalue.filter, blast.results.m.evalue.filter[,1])
  
  print(paste("Assigning taxonomy",sep=""))
  for (otu in mixedsort(names(blast.results.m.evalue.filter.list)))
  {
    print(otu)
    blast.results.m.tmp<-NULL
    accession.genus.species.m.tmp<-NULL
    sseqid.count<-NULL
    
    blast.results.m.tmp <- blast.results.m.evalue.filter.list[[otu]]
    #accession.genus.species.m.tmp <- matrix(unlist(strsplit(as.vector(blast.results.m.tmp[,2]),"_")),ncol=3,byrow=TRUE)
    accession.genus.species.m.tmp <- as.data.frame(apply(as.matrix(blast.results.m.tmp),1, function(x) unlist(strsplit(x[2],"_"))[1]))
    accession.genus.species.m.tmp <- cbind(accession.genus.species.m.tmp,as.data.frame(apply(as.matrix(blast.results.m.tmp),1, function(x) unlist(strsplit(x[2],"_"))[2])))
    accession.genus.species.m.tmp <- cbind(accession.genus.species.m.tmp,as.data.frame(apply(as.matrix(blast.results.m.tmp),1, function(x) unlist(strsplit(x[2],"_"))[3])))
    accession.genus.species.m.tmp <- cbind(accession.genus.species.m.tmp,paste(accession.genus.species.m.tmp[,2],accession.genus.species.m.tmp[,3],sep="_"))
  
    colnames(accession.genus.species.m.tmp) <- c("accession","genus","species","genusspecies")
    sseqid.count <- ddply(accession.genus.species.m.tmp,.(genusspecies),summarize,cfreq=length(genusspecies))
    print(head(sseqid.count))
    sseqid.count <- sseqid.count[order(sseqid.count$cfreq, decreasing = T),]
    sseqid.count[,2] <- sseqid.count[,2]/dim(blast.results.m.tmp)[1]

    if((sseqid.count[1,2]==1 && sseqid.count[1,1]=="Unknown_Unknown") || (sseqid.count[1,2]>=HIT_ABUNDANCE_CUTOFF && sseqid.count[1,1]=="Unknown_Unknown")){
      if(is.null(taxanomic.assignment.m)){
        tmp.m<-NULL
        tmp.m<-c("Unknown","Unknown","Unknown","Unknown","Unknown","Unknown","Unknown","Unknown")
        dim(tmp.m)<-c(1,8)
        taxanomic.assignment.m <- cbind(otu,"Unknown",tmp.m)
        colnames(taxanomic.assignment.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
      }else{
        tmp.m<-NULL
        tmp.m<-c("Unknown","Unknown","Unknown","Unknown","Unknown","Unknown","Unknown","Unknown")
        dim(tmp.m)<-c(1,8)
        tmp.m <- cbind(otu,"Unknown",tmp.m)
        colnames(tmp.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
        taxanomic.assignment.m <- rbind(taxanomic.assignment.m,tmp.m)
      } 
    }else if(sseqid.count[1,2]==1){
      if(is.null(taxanomic.assignment.m)){
        taxanomic.assignment.m <- cbind(otu,gsub(".*_","",sseqid.count[1,1]),taxanomic.csv.m[which(taxanomic.csv.m[,1]==gsub("_.*","",sseqid.count[1,1])),])
        colnames(taxanomic.assignment.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
      }else{
        tmp.m<-NULL
        tmp.m <- cbind(otu,gsub(".*_","",sseqid.count[1,1]),taxanomic.csv.m[which(taxanomic.csv.m[,1]==gsub("_.*","",sseqid.count[1,1])),])
        colnames(tmp.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
        taxanomic.assignment.m <- rbind(taxanomic.assignment.m,tmp.m)
      }
      
    }else if(sseqid.count[1,2]>=HIT_ABUNDANCE_CUTOFF){
      if(is.null(taxanomic.assignment.m)){
        taxanomic.assignment.m <- cbind(otu,gsub(".*_","",sseqid.count[1,1]),taxanomic.csv.m[which(taxanomic.csv.m[,1]==gsub("_.*","",sseqid.count[1,1])),])
        colnames(taxanomic.assignment.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
      }else{
        tmp.m<-NULL
        tmp.m <- cbind(otu,gsub(".*_","",sseqid.count[1,1]),taxanomic.csv.m[which(taxanomic.csv.m[,1]==gsub("_.*","",sseqid.count[1,1])),])
        colnames(tmp.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
        taxanomic.assignment.m <- rbind(taxanomic.assignment.m,tmp.m)
      }
      
    }else{
      all.taxa.lineage<-NULL
      all.taxa.lineage <- taxanomic.csv.m[which(taxanomic.csv.m[,1] %in% unique(gsub("_.*","",sseqid.count[,1]))),]
      if(is.null(taxanomic.assignment.m)){
        tmp.m<-NULL
        tmp.m<-apply(all.taxa.lineage,2, function(x){ifelse(length(as.vector(unique(x)))==1,unique(x),"Ambiguous")})
        dim(tmp.m)<-c(1,8)
        taxanomic.assignment.m <- cbind(otu,"Ambiguous",tmp.m)
        colnames(taxanomic.assignment.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
      }else{
        tmp.m<-NULL
        tmp.m<-apply(all.taxa.lineage,2, function(x){ifelse(length(as.vector(unique(x)))==1,unique(x),"Ambiguous")})
        dim(tmp.m)<-c(1,8)
        tmp2.m<-NULL
        tmp2.m<-cbind(otu,"Ambiguous",tmp.m)
        colnames(tmp2.m)<-c("otu","Species","Genus","Family","Order","Subclass","Class","Subphylum","Phylum", "Kingdom")
        taxanomic.assignment.m <- rbind(taxanomic.assignment.m,tmp2.m)
      }
    }
  }
  taxanomic.assignment.m <- apply(taxanomic.assignment.m,2, function(x){gsub("Incertae sedis","Unknown",x)})
  taxanomic.assignment.m <- as.data.frame(as.matrix(taxanomic.assignment.m))
  

  print(paste("Creating OTU table",sep=""))
  ucf.results.m<-as.data.frame(read.delim(file=paste(FILEDIR, ucf, sep="/"),
                                          header=F,as.is=T,row.names=NULL,fill=F,sep="\t",
                                          quote="",comment.char="",blank.lines.skip=F,stringsAsFactors=F),stringsAsFactors=F)
  
  colnames(ucf.results.m) <- c("record_type","cluster_number","sequence_length","percent_identity","strand","blank_1","blank_2","compressed_alignment","query_id","subject_id")
  
  
  ucf.results.subset.m <- ucf.results.m[which(ucf.results.m[,"record_type"]=="H"),c("record_type","cluster_number","sequence_length","percent_identity","strand","compressed_alignment","query_id","subject_id")]
  
  
  ucf.results.subset.m<-cbind(ucf.results.subset.m,gsub(".*=","",gsub(":.*","",ucf.results.subset.m[,"query_id"])))
  ucf.results.subset.m<-cbind(ucf.results.subset.m,gsub(";","",gsub(".*;size=","",ucf.results.subset.m[,"query_id"])))
  ucf.results.subset.m<-cbind(ucf.results.subset.m,gsub(";.*","",ucf.results.subset.m[,"subject_id"]))
  ucf.results.subset.m<-cbind(ucf.results.subset.m,gsub(";","",gsub(".*;size=","",ucf.results.subset.m[,"subject_id"])))
  
  colnames(ucf.results.subset.m)[9:12]<-c("sample_id","sample_raw_read_count","OTU_ID","raw_read_count")
  
  ucf.results.subset.m[as.vector(grep("barcodelabel=.*",ucf.results.subset.m[,"sample_raw_read_count"])),"sample_raw_read_count"] <- NA
  ucf.results.subset.m$"sample_raw_read_count"<-as.numeric(ucf.results.subset.m$"sample_raw_read_count")
  ucf.results.subset.m[as.vector(is.na(ucf.results.subset.m[,"sample_raw_read_count"])),"sample_raw_read_count"] <- 1
  
  print(paste("Getting raw counts per sample for each OTU",sep=""))
  ucf.results.subset.otu.sample.counts <- ddply(as.data.frame(ucf.results.subset.m),c("sample_id","OTU_ID"),summarise,sum=sum(as.numeric(sample_raw_read_count)))
  
  ucf.results.subset.otu.total.counts <- ddply(ucf.results.subset.m,.(OTU_ID),summarize,total_count=unique(raw_read_count))
  
  ucf.results.subset.otu.m <- as.data.frame(matrix(data=0,nrow=length(as.vector(unique(ucf.results.subset.m[,"OTU_ID"]))),
                                                   ncol=length(as.vector(unique(ucf.results.subset.m[,"sample_id"])))),
                                            row.names=mixedsort(as.vector(unique(ucf.results.subset.m[,"OTU_ID"]))))
  colnames(ucf.results.subset.otu.m) <- mixedsort(as.vector(unique(ucf.results.subset.m[,"sample_id"])))
  ucf.results.subset.otu.sample.counts.list <- split(ucf.results.subset.otu.sample.counts[,2:3], ucf.results.subset.otu.sample.counts[,1])
  
  for(sampleid in mixedsort(names(ucf.results.subset.otu.sample.counts.list))){
    otuid.tmp <- ucf.results.subset.otu.sample.counts.list[[sampleid]]
    ucf.results.subset.otu.m[as.vector(otuid.tmp[,"OTU_ID"]),sampleid]<-otuid.tmp[,"sum"]
  }
  
  print(paste("Assigning taxanomy for each OTU",sep=""))
  rownames(taxanomic.assignment.m) <- gsub(";.*","",taxanomic.assignment.m[,1])
  
  ucf.results.subset.otu.taxa.m <- cbind(gsub(";.*","",taxanomic.assignment.m[,1]),paste("k__",taxanomic.assignment.m[,10],"; ","p__",taxanomic.assignment.m[,9],"; ",
                                                                                         "c__",taxanomic.assignment.m[,7],"; ","o__",taxanomic.assignment.m[,5],"; ",
                                                                                         "f__",taxanomic.assignment.m[,4],"; ","g__",taxanomic.assignment.m[,3],"; ",
                                                                                         "s__",taxanomic.assignment.m[,2],"; ",sep=""))
  
  unassigned.tmp <- cbind(rownames(ucf.results.subset.otu.m)[-(which(rownames(ucf.results.subset.otu.m) %in% rownames(taxanomic.assignment.m)))],paste("k__Unassigned; ","p__Unassigned; ",
                                                                                                                                                       "c__Unassigned; ","o__Unassigned; ",
                                                                                                                                                       "f__Unassigned; ","g__Unassigned; ",
                                                                                                                                                       "s__Unassigned; ",sep=""))
  unassigned.taxa.m <- as.data.frame(matrix(data="Unassigned",nrow=dim(unassigned.tmp)[1],ncol=dim(taxanomic.assignment.m)[2],dimnames=list(unassigned.tmp[,1],colnames(taxanomic.assignment.m))))
  unassigned.taxa.m[,1] <- unassigned.tmp[,1]
  
  taxanomic.assignment.m <- rbind(taxanomic.assignment.m,unassigned.taxa.m)
  
  
  names(taxanomic.assignment.m$otu) <- "NULL"
  names(taxanomic.assignment.m$Species) <- "NULL"
  names(taxanomic.assignment.m$Genus) <- "NULL"
  names(taxanomic.assignment.m$Family) <- "NULL"
  names(taxanomic.assignment.m$Order) <- "NULL"
  names(taxanomic.assignment.m$Subclass) <- "NULL"
  names(taxanomic.assignment.m$Class) <- "NULL"
  names(taxanomic.assignment.m$Subphylum) <- "NULL"
  names(taxanomic.assignment.m$Phylum) <- "NULL"
  names(taxanomic.assignment.m$Kingdom) <- "NULL"
  
  print(paste("Getting taxonomic counts",sep=""))
  species.count.m <- ddply(taxanomic.assignment.m,.(Species),"nrow")
  genus.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Genus),"nrow")
  family.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Family),"nrow")
  order.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Order),"nrow")
  subclass.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Subclass),"nrow")
  class.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Class),"nrow")
  subphylum.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Subphylum),"nrow")
  phylum.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Phylum),"nrow")
  kingdom.count.m <- ddply(data.frame(taxanomic.assignment.m),.(Kingdom),"nrow")
  
  print(paste("Creating taxonomic files",sep=""))
  dir.create(file.path(FILEDIR, "output"), showWarnings = FALSE)
  
  write.table(taxanomic.assignment.m,file=paste(FILEDIR, "output", gsub("\\..*",".results",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  
  write.table(species.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_species.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(genus.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_genus.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(family.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_family.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(order.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_order.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(subclass.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_subclass.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(class.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_class.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(subphylum.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_subphylum.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(phylum.count.m,file=paste(FILEDIR, "output", gsub("\\..*","_phylum.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  write.table(kingdom.count.m,file=paste(FILEDIR, "output", sub("\\..*","_kingdom.txt",bf,perl = T),sep="/"),row.names=F,quote=F, sep = "\t")
  


  
  ucf.results.subset.otu.taxa.m <- rbind(ucf.results.subset.otu.taxa.m,unassigned.tmp)
  
  print(paste("Creating OTU taxanomy file",sep=""))
  write.table(ucf.results.subset.otu.taxa.m,file=paste(FILEDIR, "output", sub("\\..*","_biom_taxa.txt",bf,perl = T),sep="/"),row.names=F,
              col.names=c("#OTUID","taxonomy"),quote=F, sep = "\t")
  
  print(paste("Creating BIOM text file",sep=""))
  write.table(cbind(rownames(ucf.results.subset.otu.m),ucf.results.subset.otu.m),file=paste(FILEDIR, "output", sub("\\..*","_biom.txt",bf,perl = T),sep="/"),
              row.names=F,col.names=c("#OTU ID",colnames(ucf.results.subset.otu.m)),quote=F, sep = "\t")
}