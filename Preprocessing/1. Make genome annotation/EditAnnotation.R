###########################################Part 1########################################

library("ape"); library(dplyr); library("DescTools")

options(scipen=999)

setwd("/media/claired/906143b9-524c-4abe-8428-6fc9300eb700/LatestAnnotation_June25/")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#1. work out multiple transcript "issue"
gtfHs <- rtracklayer::import('/media/claired/906143b9-524c-4abe-8428-6fc9300eb700/LatestAnnotation_June25/Homo_sapiens.GRCh38.114.gtf') #get your genome annotation, this code has specifically been deigned for Human annotations please adjust accordingly for other species
gtf_HsDf <- as.data.frame(gtfHs)
dim(gtf_HsDf[gtf_HsDf$type == "exon",])
dim(gtf_HsDf[gtf_HsDf$type == "gene",])
dim(gtf_HsDf[gtf_HsDf$type == "transcript",])

unique(levels(gtf_HsDf$type)) 
levels(gtf_HsDf$seqnames) 
unique(gtf_HsDf$gene_biotype)
unique(gtf_HsDf$transcript_biotype)
length(unique(gtf_HsDf[gtf_HsDf$type == "exon",]$transcript_id))

#get genes with multiple  transcript "child" features
length(unique(gtf_HsDf[gtf_HsDf$type == "transcript",]$gene_id))
#get gene IDs with multiple transcripts
ParentInfo <- data.frame(gtf_HsDf[gtf_HsDf$type == "transcript",]$gene_id, gtf_HsDf[gtf_HsDf$type == "transcript",]$transcript_id)
colnames(ParentInfo) <- c("Gene_ID", "Transcript_ID")
length(unique(ParentInfo$Gene_ID))
length(unique(ParentInfo$Transcript_ID))

occurances <- c()
for (parent in unique(ParentInfo$Gene_ID)) {
  occurances <- append(occurances, nrow(ParentInfo[ParentInfo$Gene_ID == parent,]))
}
parentinfo <- data.frame(unique(ParentInfo$Gene_ID), occurances)
colnames(parentinfo) <- c("Parent", "Children")
dim(parentinfo[parentinfo[,2] == 1,]) #44,630
dim(parentinfo[parentinfo[,2] != 1,]) #34,264
dim(parentinfo[parentinfo[,2] == 2,]) #6,674
dim(parentinfo[parentinfo[,2] == 3,]) #4,239
dim(parentinfo[parentinfo[,2] > 3,]) #23,351
MultipleTranscripts <- parentinfo[parentinfo[,2] != 1,1]

write.csv(MultipleTranscripts, "MultipleTRanscriptGene.csv")

#use gene_id to collect exons, exon_number will need to be altered, ignore CDS
exons <- gtf_HsDf[gtf_HsDf$type == "exon",]
CDS <- gtf_HsDf[gtf_HsDf$type == "CDS",]
New_gtf <- data.frame(matrix(ncol=27,nrow=0), stringsAsFactors = FALSE)
ToDelExon <- c()
ToDELTranscripts <- c()
AllTransDontOverlap <- c()
#touching_exons <- c()
SaveInfoUsed <- data.frame(matrix(ncol=5,nrow=0), stringsAsFactors = FALSE) #Gene, MegaName, number of transcripts used, TranscriptIDs, Transcript Biotypes
for (parent in MultipleTranscripts) { #for each gene with multiple transcripts...
  #get transcripts
  children <- gtf_HsDf[gtf_HsDf$type == "transcript" & gtf_HsDf$gene_id == parent,]
  originalChildrenNum <- nrow(children)
  #only get transcripts that overlap each other
  children_red <- children %>% 
    arrange(start) %>% 
    group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
    summarise(start = first(start), end = max(end))
  NewChildrenNum <- nrow(children_red)
  #depending on the number of new relative to old children:
  #same number - therefore none overlap - count these so you know how many you have later but ultimately ignore
  if (NewChildrenNum == originalChildrenNum) {
    AllTransDontOverlap <- append(AllTransDontOverlap, parent)
  }
  #1 - take all forward
  else if (NewChildrenNum == 1) {
    width <- children_red$end-children_red$start+1
    myNewtranscriptID <- paste(parent, "_Mega", sep="")
    #make new transcript row
    TranscriptRowToEdit <- children[1,]
    TranscriptRowToEdit$start <- children_red$start
    TranscriptRowToEdit$end <- children_red$end
    TranscriptRowToEdit$width <- width
    TranscriptRowToEdit$transcript_id <- myNewtranscriptID
    info <- c(parent, myNewtranscriptID, nrow(children), paste(unique(children$transcript_id), collapse="_"), paste(unique(children$transcript_biotype), collapse="_"))
    SaveInfoUsed <- rbind(SaveInfoUsed, info, stringsAsFactors=FALSE)
    colnames(SaveInfoUsed) <- c("gene_id", "NewTranscriptID", "NumberofTranscripts", "OriginalTranscriptIds", "OriginalTranscriptBiotype")
    #get exons
    MyExons <- exons[exons$transcript_id %in% children$transcript_id,]
    MyExonsRed <- MyExons %>% 
      arrange(start) %>% 
      group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
      summarise(start = first(start), end = max(end))
    #MyExonsRed_test <- MyExonsRed
    #MyExonsRed_test[,2] <- MyExonsRed_test[,2] - 1
    #MyExonsRed_test[,3] <- MyExonsRed_test[,3] + 1
    #  MyExonsRed_Further <- MyExonsRed_test %>% 
    #  arrange(start) %>% 
    #  group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
    #  summarise(start = first(start), end = max(end))
    #if (nrow(MyExonsRed_Further) < nrow(MyExonsRed_test)) {
    #  touching_exons <- append(touching_exons, parent)
    #find overlap with MyExonsRed_test and take outtermost original coordinates!
    #}
    ID <- c()
    for (exon in 1:nrow(MyExonsRed)) {
      ID <- append(ID, paste(myNewtranscriptID, "-E", exon, sep=""))
    }
    ExNumber <- c(1:nrow(MyExonsRed))
    start <- MyExonsRed[,2]
    end <- MyExonsRed[,3]
    width <- end - start + 1
    start <- start %>% pull (start)
    end <- end %>% pull (end)
    width <- width %>% pull (end)
    if (MyExons$strand[1] == "-") {
      start <- rev(start)
      end <- rev(end)
      width <- rev(width)
    }
    #make new exon rows 
    ExonRowstoEdit <- MyExons[1:nrow(MyExonsRed),]
    ExonRowstoEdit$start <- start
    ExonRowstoEdit$end <- end
    ExonRowstoEdit$width <- width
    ExonRowstoEdit$transcript_support_level <- TranscriptRowToEdit$transcript_support_level
    ExonRowstoEdit$tag <- TranscriptRowToEdit$tag
    ExonRowstoEdit$transcript_version <- TranscriptRowToEdit$transcript_version
    ExonRowstoEdit$transcript_name <- TranscriptRowToEdit$transcript_name
    ExonRowstoEdit$transcript_source <- TranscriptRowToEdit$transcript_source
    ExonRowstoEdit$transcript_biotype <- TranscriptRowToEdit$transcript_biotype
    ExonRowstoEdit$ccds_id <- TranscriptRowToEdit$ccds_id
    ExonRowstoEdit$transcript_id <- rep(myNewtranscriptID, nrow(MyExonsRed))
    ExonRowstoEdit$exon_number <- ExNumber
    ExonRowstoEdit$exon_id <- ID
    #save together
    NewRows <- rbind(TranscriptRowToEdit, ExonRowstoEdit)
    New_gtf <- rbind(New_gtf, NewRows)
    #lines to remove
    ToDelExon <- append(ToDelExon, MyExons$exon_id)
    ToDELTranscripts <- append(ToDELTranscripts, children$transcript_id)
  }
  #less but not equal to 1 - take those that overlap forward, 
  else if (NewChildrenNum > 1 & NewChildrenNum < originalChildrenNum) {
    for (newChild in 1:nrow(children_red)) {
      mychild <- children_red[newChild,]
      width <- mychild$end-mychild$start+1
      myNewtranscriptID <- paste(parent, "_", toString(newChild), "_Mega", sep="")
      #get the transcripts that fit within the new ones and use those details for next steps
      overlapping <- children[children$start >= mychild$start & children$end <= mychild$end,] #due to nature they have to be within?
      myTranscripts <- overlapping$transcript_id
      #just take transcript info from first occurance
      #make new transcript row
      ImportantTranscriptRows <- children[children$transcript_id %in% myTranscripts, ]
      TranscriptRowToEdit <- ImportantTranscriptRows[1,]
      TranscriptRowToEdit$start <- mychild$start
      TranscriptRowToEdit$end <- mychild$end
      TranscriptRowToEdit$width <- width
      TranscriptRowToEdit$transcript_id <- myNewtranscriptID
      info <- c(parent, myNewtranscriptID, nrow(ImportantTranscriptRows), paste(unique(myTranscripts), collapse="_"), paste(unique(ImportantTranscriptRows$transcript_biotype), collapse="_"))
      SaveInfoUsed <- rbind(SaveInfoUsed, info, stringsAsFactors=FALSE)
      colnames(SaveInfoUsed) <- c("gene_id", "NewTranscriptID", "NumberofTranscripts", "OriginalTranscriptIds", "OriginalTranscriptBiotype")
      #get exons
      MyExons <- exons[exons$transcript_id %in% myTranscripts,]
      MyExonsRed <- MyExons %>% 
        arrange(start) %>% 
        group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
        summarise(start = first(start), end = max(end))
      #MyExonsRed_test <- MyExonsRed
      #MyExonsRed_test[,2] <- MyExonsRed_test[,2] - 1
      #MyExonsRed_test[,3] <- MyExonsRed_test[,3] + 1
      #MyExonsRed_Further <- MyExonsRed_test %>% 
      #  arrange(start) %>%  example:ENSG00000205670
      #  group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
      #  summarise(start = first(start), end = max(end))
      # if (nrow(MyExonsRed_Further) < nrow(MyExonsRed_test)) {
      #   touching_exons <- append(touching_exons, parent)
      #find overlap with MyExonsRed_test and take outtermost original coordinates!
      #}
      ID <- c()
      for (exon in 1:nrow(MyExonsRed)) {
        ID <- append(ID, paste(myNewtranscriptID, "-E", exon, sep=""))
      }
      ExNumber <- c(1:nrow(MyExonsRed))
      start <- MyExonsRed[,2]
      end <- MyExonsRed[,3]
      width <- end - start + 1
      start <- start %>% pull (start)
      end <- end %>% pull (end)
      width <- width %>% pull (end)
      if (MyExons$strand[1] == "-") {
        start <- rev(start)
        end <- rev(end)
        width <- rev(width)
      }
      #make new exon rows 
      ExonRowstoEdit <- MyExons[1:nrow(MyExonsRed),]
      ExonRowstoEdit$start <- start
      ExonRowstoEdit$end <- end
      ExonRowstoEdit$width <- width
      ExonRowstoEdit$transcript_support_level <- TranscriptRowToEdit$transcript_support_level
      ExonRowstoEdit$tag <- TranscriptRowToEdit$tag
      ExonRowstoEdit$transcript_version <- TranscriptRowToEdit$transcript_version
      ExonRowstoEdit$transcript_name <- TranscriptRowToEdit$transcript_name
      ExonRowstoEdit$transcript_source <- TranscriptRowToEdit$transcript_source
      ExonRowstoEdit$transcript_biotype <- TranscriptRowToEdit$transcript_biotype
      ExonRowstoEdit$ccds_id <- TranscriptRowToEdit$ccds_id
      ExonRowstoEdit$transcript_id <- rep(myNewtranscriptID, nrow(MyExonsRed))
      ExonRowstoEdit$exon_number <- ExNumber
      ExonRowstoEdit$exon_id <- ID
      #save together
      NewRows <- rbind(TranscriptRowToEdit, ExonRowstoEdit)
      New_gtf <- rbind(New_gtf, NewRows)
      #lines to remove
      ToDelExon <- append(ToDelExon, MyExons$exon_id)
      ToDELTranscripts <- append(ToDELTranscripts, myTranscripts)
    }
  }
}
options(scipen=999)
#print(count)
print(AllTransDontOverlap) #ENSG00000272027
print(length(AllTransDontOverlap)) #1
#print(length(touching_exons)) #0
write.csv(SaveInfoUsed, "MultipleTRanscriptGeneDetails.csv", row.names=FALSE)
write.csv(New_gtf, "MegaTranscriptsGTFlines.csv", row.names=FALSE)
write.csv(AllTransDontOverlap, "MultipleTranscriptButNoOverlap.csv", row.names=FALSE)
write.csv(ToDELTranscripts, "ToDELTranscripts.csv", row.names=FALSE)

#New_gtf <- read.csv("MegaTranscriptsGTFlines.csv")
#AllTransDontOverlap <- read.csv("MultipleTranscriptButNoOverlap.csv")

#non-overlapping at all transcripts
length(unique(gtf_HsDf[gtf_HsDf$type == "transcript",]$transcript_id)) #38,7954
length(unique(New_gtf[New_gtf$type == "transcript",]$transcript_id)) #34,263 -1 less than previous estimate
length(unique(gtf_HsDf[gtf_HsDf$type == "gene",]$gene_id)) #78,894
#unequal number of parents and genes is due to genes with multiple non-overlapping genes, if these were non-overlapping in the original file they keep their original IDs, if they were reduced down from larger sets then they get new IDs with "_transcript number"
New_gtf[New_gtf$type == "transcript" & grepl("_1_Mega", New_gtf$transcript_id),] #none
New_gtf[New_gtf$type == "transcript" & grepl("_2_Mega", New_gtf$transcript_id),] #none

'%ni%' <- Negate('%in%')
MultipleTranscripts <- MultipleTranscripts[MultipleTranscripts %ni% AllTransDontOverlap] #

#check sum of these + number of transcripts match the number of lines removed at next step
nrow(gtf_HsDf[gtf_HsDf$type == "three_prime_utr" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +
  nrow(gtf_HsDf[gtf_HsDf$type == "five_prime_utr" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +
  nrow(gtf_HsDf[gtf_HsDf$type == "start_codon" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +#
  nrow(gtf_HsDf[gtf_HsDf$type == "stop_codon" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +#
  nrow(gtf_HsDf[gtf_HsDf$type == "CDS" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +#
  nrow(gtf_HsDf[gtf_HsDf$type == "exon" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +#
  nrow(gtf_HsDf[gtf_HsDf$type == "Selenocysteine" & gtf_HsDf$transcript_id %in% ToDELTranscripts,]) +#
  length(ToDELTranscripts)
length(unique(ToDELTranscripts)) #sanity check that these are the same

gtf_HsDfEd <- gtf_HsDf[gtf_HsDf$transcript_id %ni% ToDELTranscripts,] #
nrow(gtf_HsDf)-nrow(gtf_HsDfEd)
#remove CDS
gtf_HsDfEd <- gtf_HsDfEd[gtf_HsDfEd$type != "CDS",]

#put in to df with ee after gene lines
TopLevelParents <- MultipleTranscripts
existingDF <- gtf_HsDfEd
rownames(existingDF) <- 1:nrow(existingDF)
nrow(existingDF) + nrow(New_gtf) #614161
for (Parent in unique(TopLevelParents)) {
  parent_row <- existingDF[existingDF$gene_id == Parent,] #get the genes row
  parent_row <- parent_row[complete.cases(parent_row[,1:3]),]
  position <- rownames(parent_row) #get its position
  new_rows <- New_gtf[grep(Parent, New_gtf$gene_id), ] #get rows with this parent from the mega transcript info: transcripts and exons
  if (nrow(new_rows) > 0) {
    existingDF <- rbind(existingDF[1:position,], new_rows, existingDF[-(1:position),])
    rownames(existingDF) <- 1:nrow(existingDF)   
  }
}
dim(existingDF) #check its old + nrow of New_gtf

options(scipen=999)
write.csv(existingDF, "MegaAddedToNormal.csv", row.names=FALSE)

unique(existingDF$type)
existingDF[existingDF$type == "Selenocysteine",] #2 nt bits, 119 of them, are tied to genes and transcripts
#remove these
existingDF <- existingDF[existingDF$type != "Selenocysteine",]

#check number of unique parents is genes + multi-transcript etc.
#Get the number/proportion of single exon genes
existingDFExons <-  existingDF[existingDF$type == "exon",] #
ExPa <- data.frame(existingDFExons$exon_id, existingDFExons$transcript_id)
colnames(ExPa) <- c("ID", "Parent")
length(unique(ExPa$Parent))
occurancesExPa <- c()
for (parent in unique(ExPa$Parent)) {
  occurancesExPa <- append(occurancesExPa, nrow(ExPa[ExPa$Parent == parent,]))
}
parentinfoExPa <- cbind(unique(ExPa$Parent), occurancesExPa)
colnames(parentinfoExPa) <- c("Parent", "Children")
nrow(parentinfoExPa[parentinfoExPa[,2] == 1,])/nrow(parentinfoExPa)*100
nrow(parentinfoExPa[parentinfoExPa[,2] >= 2,])/nrow(parentinfoExPa)*100
nrow(parentinfoExPa[parentinfoExPa[,2] >= 3,])/nrow(parentinfoExPa)*100 

#change exon IDs for exons for other transcripts - this line should laways be checked to ensure the order of columns has not changed:
existingDF[existingDF$type == "exon" & !grepl("Mega", existingDF$transcript_id), 23] <- paste(existingDF[existingDF$type == "exon" & !grepl("Mega", existingDF$transcript_id), 15], "-E", existingDF[existingDF$type == "exon" & !grepl("Mega", existingDF$transcript_id), 22], sep="")

options(scipen=999)
write.csv(existingDF, "All_withFixedExonID.csv", row.names=FALSE)

#check number of parents
length(unique(existingDF[!grepl("Mega", existingDF$transcript_id), 15]))+
  length(unique(existingDF[grepl("Mega", existingDF$transcript_id), 15]))
length(unique(existingDF[grepl("_1_Mega", existingDF$transcript_id), 15])) 
length(unique(existingDF[grepl("_2_Mega", existingDF$transcript_id), 15]))

setdiff(c(unique(existingDF[!grepl("Mega", existingDF$transcript_id), 15]),unique(existingDF[grepl("Mega", existingDF$transcript_id), 15])), unique(ExPa$Parent))

#save as gtf
colnames(existingDF) #add dots and reorder etc. #here needs to be written onwards
NewGffFile <- existingDF
NewGffFile <- NewGffFile[,-4] #remove width column
dim(NewGffFile) #614160 26
#seqnames source type start end . strand . merged:gene_id gene_name gene_biotype transcript
NewGffFileBeg <- NewGffFile[,c(1,5,6,2,3,4)]
NewGffFileBeg$dot1 <- rep(".", nrow(NewGffFileBeg))
NewGffFileBeg$dot2 <- rep(".", nrow(NewGffFileBeg))
NewGffFileBeg <- NewGffFileBeg[,c(1,2,3,4,5,7,6,8)]
NewGffFileBeg2 <- NewGffFileBeg
#sort out the final merge columm, remove NA, merge exon_id gene_id transcript_id, parent
finalset <- NewGffFile[, c(9,10,11,12, 13,14,15,16,17,18,19,20,21,22,23,24,25,26)]
finalset2 <- NewGffFile[, c(6,9,10,11,12, 13,14,15,16,17,18,19,20,21,22,23,24,25,26)]

#Sort out ID information for features
ID <- c()
Parent <- c()
for (line in 1:nrow(finalset2)){
  MyLine <- finalset2[line,]
  if (MyLine$type == "gene") {
    ID <- append(ID,MyLine[,2])
    Parent <- append(Parent,"NA")
  }
  else if (MyLine$type == "transcript") {
    ID <- append(ID,MyLine[,7])
    Parent <- append(Parent,MyLine[,2])
  }
  else if (MyLine$type == "exon") {
    ID <- append(ID,MyLine[,15]) #other line changed to fix exon ID issue
    Parent <- append(Parent,MyLine[,7])
  }
  else if (MyLine$type == "start_codon") {
    ID <- append(ID,paste("start_codon_", MyLine[,7], sep=""))
    Parent <- append(Parent,MyLine[,7])
  }
  else if (MyLine$type == "stop_codon") {
    ID <- append(ID,paste("stop_codon_", MyLine[,7], sep=""))
    Parent <- append(Parent,MyLine[,7])
  }
  else if (MyLine$type == "five_prime_utr") {
    ID <- append(ID,paste("five_prime_utr_", MyLine[,7], sep=""))
    Parent <- append(Parent,MyLine[,7])
  }
  else if (MyLine$type == "three_prime_utr") {
    ID <- append(ID,paste("three_prime_utr_", MyLine[,7], sep=""))
    Parent <- append(Parent,MyLine[,7])
  }
}
finalset3 <- data.frame(ID, Parent, finalset2)
finalset3<- finalset3[,-3]
for (column_name in 1:ncol(finalset3)) { #combine: ID, DEScription, Ontology_term, Parent, protein_source_id, etc. with ; between them into new final row: callled everythignelse
  name <- paste(colnames(finalset3)[column_name], "=", sep="")
  finalset3[, column_name] <- glue::glue("{name}{finalset3[,column_name]}")
}
for (column_name in 1:19) {
  finalset3[,column_name] <- paste0(finalset3[,column_name], ";")
}
merged <- do.call(paste0, finalset3[colnames(finalset3)])
NewGffFileBeg$merged <- merged
write.csv(NewGffFileBeg, 'Homo_sapiens.GRCh38.104_EditedGff3.csv', row.names=FALSE)

#for gtf
for (column_name in 1:ncol(finalset)) { #combine: ID, DEScription, Ontology_term, Parent, protein_source_id, etc. with ; between them into new final row: callled everythignelse
  finalset[,column_name] <- paste('"', finalset[,column_name], '"', sep="")
  name <- paste0(colnames(finalset)[column_name], " ", sep="")
  finalset[, column_name] <- glue::glue("{name}{finalset[,column_name]}")
}
for (column_name in 1:17) {
  finalset[,column_name] <- paste0(finalset[,column_name], "; ")
}
merged2 <- do.call(paste0, finalset[colnames(finalset)])
write.csv(data.frame(NewGffFileBeg2, merged2), 'Homo_sapiens.GRCh38.104_EditedGTF.csv', row.names=FALSE)
write.table(data.frame(NewGffFileBeg2, merged2), 'Homo_sapiens.GRCh38.104_Edited.gtf', row.names=FALSE, sep = '\t',  col.names=FALSE, quote = FALSE)
#Fix the Parents:
FixNA <- read.csv("Homo_sapiens.GRCh38.104_EditedGff3.csv")
FixNA$merged <- gsub("Parent=NA;", "", FixNA$merged)
FixNA$source <- rep("ensembl_havana", 614160) #check this number for future annoatation reference versions
write.table(FixNA, 'Homo_sapiens.GRCh38.104_Edited.gff3', row.names=FALSE, sep = "\t",  col.names=FALSE, quote = FALSE)

#########################Part 2#################################

#assign intron IDs
options(scipen=999) #make sure they are in the inverse order for "-" stranded genes
gffIn <- rtracklayer::import('Homo_sapiens.GRCh38.104_Edited_IntronsAdded.gff3')
gffIn_df <- as.data.frame(gffIn)
dim(unique(gffIn_df[gffIn_df$type == "gene", ]))
dim(unique(gffIn_df[gffIn_df$type == "exon", ]))
dim(unique(gffIn_df[gffIn_df$type == "transcript", ]))
dim(unique(gffIn_df[gffIn_df$type == "intron", ]))

#examine length distributions
plot(density(gffIn_df[gffIn_df$type == "exon",]$width))
plot(density(gffIn_df[gffIn_df$type == "intron",]$width))

intron_parents <- unique(gffIn_df[gffIn_df$type == "intron",]$Parent) #
introns <- gffIn_df[gffIn_df$type == "intron",] #
exons <- gffIn_df[gffIn_df$type == "exon",] #

for (parent in intron_parents) { #assign intron IDs
  Myrows <- introns[introns$Parent == parent,]
  MyrowsEx <- exons[exons$Parent == parent,]
  numberofintrons <- nrow(Myrows)
  IntronIDs <- c()
  for (intronNum in 1:numberofintrons) {
    IntronIDs <- append(IntronIDs, paste(parent, "-I", toString(intronNum), sep=""))
  }
  if (Myrows$strand[1] == "-") {
    IntronIDs <- rev(IntronIDs)
  }
  Myrows[, 10] <- IntronIDs
  Myrows[, 6] <- rep(MyrowsEx[1, 6], numberofintrons)
  Myrows[, 11] <- rep(MyrowsEx[1, 11], numberofintrons)
  Myrows[, 16] <- rep(MyrowsEx[1, 16], numberofintrons)
  Myrows[, 12] <- rep(MyrowsEx[1, 12], numberofintrons)
  Myrows[, 13] <- rep(MyrowsEx[1, 13], numberofintrons)
  Myrows[, 14] <- rep(MyrowsEx[1, 14], numberofintrons)
  Myrows[, 15] <- rep(MyrowsEx[1, 15], numberofintrons)
  Myrows[, 17] <- rep(MyrowsEx[1, 17], numberofintrons)
  Myrows[, 18] <- rep(MyrowsEx[1, 18], numberofintrons)
  Myrows[, 19] <- rep(MyrowsEx[1, 19], numberofintrons)
  Myrows[, 20] <- rep(MyrowsEx[1, 20], numberofintrons)
  introns[introns$Parent == parent,] <- Myrows
}

options(scipen=999)
gtf_dfcopy <- gffIn_df
gtf_dfcopy[gtf_dfcopy$type == "intron",] <- introns

#reverse order for negatively stranded introns 
NegStrand <- gtf_dfcopy[gtf_dfcopy$strand == "-",] #they have been reversed here
Parents <- NegStrand$transcript_id
for (transcript in Parents) {
  exons <- NegStrand[NegStrand$transcript_id == transcript & NegStrand$type == "exon",]
  if (nrow(exons) > 1) {
    Rearranged <- arrange(NegStrand[NegStrand$transcript_id == transcript & NegStrand$type != "transcript",], -row_number())
    NegStrand[NegStrand$transcript_id == transcript & NegStrand$type != "transcript",] <- Rearranged
  }
}
gtf_dfcopy[gtf_dfcopy$strand == "-",] <- NegStrand #774269

#write to gff and gtf
colnames(gtf_dfcopy) #add dots and reorder etc. #here needs to be written on wards
NewGffFile <- gtf_dfcopy
NewGffFile <- NewGffFile[,-4] #remove width column
dim(NewGffFile) #

NewGffFileBeg <- NewGffFile[,c(1,5,6,2,3,4)]
NewGffFileBeg$dot1 <- rep(".", nrow(NewGffFileBeg))
NewGffFileBeg$dot2 <- rep(".", nrow(NewGffFileBeg))
NewGffFileBeg <- NewGffFileBeg[,c(1,2,3,4,5,7,6,8)]

#sort out the final merge columm, remove NA, merge exon_id gene_id transcript_id, parent
finalset <- NewGffFile[, c(9,10,11,12, 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)]

#GFF3:
finalset2 <- finalset
for (column_name in 1:ncol(finalset2)) { #combine: ID, DEScription, Ontology_term, Parent, protein_source_id, etc. with ; between them into new final row: callled everythignelse
  name <- paste(colnames(finalset2)[column_name], "=", sep="")
  finalset2[, column_name] <- glue::glue("{name}{finalset2[,column_name]}")
}
for (column_name in 1:19) {
  finalset2[,column_name] <- paste0(finalset2[,column_name], ";")
}
merged <- do.call(paste0, finalset2[colnames(finalset2)])
merged <- gsub(";Parent=character\\(0\\)", "", merged)
NewGffFileBeg2 <- NewGffFileBeg
NewGffFileBeg2$merged <- merged
write.csv(NewGffFileBeg2, 'Homo_sapiens.GRCh38.104_IntronsFinalGff3.csv', row.names=FALSE)
write.table(NewGffFileBeg2, 'Homo_sapiens.GRCh38.104_IntronsFinal.gff3', row.names=FALSE, sep = '\t',  col.names=FALSE, quote = FALSE)

NewGffFileBeg2[grepl("ENSG00000248333", NewGffFileBeg2$merged),]

#for the gtf
finalset3 <- finalset
for (column_name in 1:ncol(finalset3)) { #combine: ID, DEScription, Ontology_term, Parent, protein_source_id, etc. with ; between them into new final row: callled everythignelse
  finalset3[,column_name] <- paste('"', finalset3[,column_name], '"', sep="")
  name <- paste0(colnames(finalset3)[column_name], " ", sep="")
  finalset3[, column_name] <- glue::glue("{name}{finalset3[,column_name]}")
}
for (column_name in 1:19) {
  finalset3[,column_name] <- paste0(finalset3[,column_name], "; ")
}
merged2 <- do.call(paste0, finalset3[colnames(finalset3)])
merged2 <- gsub("; Parent \"character\\(0\\)\"", "", merged2)

write.csv(data.frame(NewGffFileBeg, merged2), 'Homo_sapiens.GRCh38.104_IntronsFinalGTF.csv', row.names=FALSE)
write.table(data.frame(NewGffFileBeg, merged2), 'Homo_sapiens.GRCh38.104_IntronsFinal.gtf', row.names=FALSE, sep = '\t',  col.names=FALSE, quote = FALSE)
