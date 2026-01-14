#Main Steps:
#1. remove non-expressed
#2. remove absent in other
#3. remove genes with multiple non-overlapping transcripts
#4. remove retained introns
#5. remove genes where introns overlap exons of other genes
#6. normalize together or separately
#7. save as output files

CalculateSplicedUnspliced <- function(Intron, Exon, Out, TranscriptList, GeneInfo, gff3file, normalise, OutName) {
  library(limma); library(edgeR); library(ggfortify); library("FactoMineR"); library("factoextra"); library("tidyr"); library("ape"); library("dplyr"); library("DescTools")
  `%ni%` = Negate(`%in%`)
  
  IntronMeta <- read.csv(Intron)
  rownames(IntronMeta) <- IntronMeta$ID
  IntronMeta <- IntronMeta[,-1]
  
  ExonMeta <- read.csv(Exon)
  rownames(ExonMeta) <- ExonMeta$ID
  ExonMeta <- ExonMeta[,-1]
  
  setwd(Out)

  #remove non-expressed
  IntronMeta <- IntronMeta[which(rowSums(IntronMeta >= 1) >= 1),] #remove non-expressed
  ExonMeta <- ExonMeta[which(rowSums(ExonMeta >= 1) >= 1),] #remove non-expressed

  #remove those absent in the other
  IntronMeta_filtered <- IntronMeta[rownames(IntronMeta) %in% rownames(ExonMeta),] #
  ExonMeta_filtered <- ExonMeta[rownames(ExonMeta) %in% rownames(IntronMeta),] #
  
  #re-order
  IntronMeta_filtered <- IntronMeta_filtered[match(rownames(ExonMeta_filtered), rownames(IntronMeta_filtered)),]
  
  print("Introns and exons are in the same order?")
  print(identical(rownames(IntronMeta_filtered), rownames(ExonMeta_filtered))) #should be true
  
  print("Samples are in the same order?")
  print(identical(colnames(IntronMeta_filtered), colnames(ExonMeta_filtered))) #should be true
  
  #remove genes with multiple non-overlapping transcripts
  IntronMeta_filtered <- IntronMeta_filtered[rownames(IntronMeta_filtered) %ni% TranscriptList,]
  ExonMeta_filtered <- ExonMeta_filtered[rownames(ExonMeta_filtered) %ni% TranscriptList,]

  print("Introns and exons are still in the same order?")
  print(identical(rownames(IntronMeta_filtered), rownames(ExonMeta_filtered))) #should be true   

  #get gene IDs - might be needed for matching with the DEGs later
  Transcript_Gene_Matches <- read.csv(GeneInfo)
  ExonGeneIDs <- gsub("_Mega", "", rownames(ExonMeta_filtered))
  IDs <- c()
  for (gene in ExonGeneIDs) {
    if (grepl("ENST", gene) == TRUE) {
      IDs <- append(IDs, as.character(Transcript_Gene_Matches[Transcript_Gene_Matches$TranscriptID == gene,]$GeneID))
    }
    else {
      IDs <- append(IDs, gene)
    }
  }
  Translation <- data.frame(ExonGeneIDs, IDs)
  print("Exons IDs are in the same order?")
  print(identical(unlist(lapply(Translation$ExonGeneIDs, toString)), unlist(lapply(ExonGeneIDs, toString))))
  write.csv(Translation, "Translation.csv", row.names=FALSE)
  
  rownames(ExonMeta_filtered) <- Translation$IDs
  rownames(IntronMeta_filtered) <- Translation$IDs

  #rm retained intron etc.
  gtfHs <- rtracklayer::import(gff3file)
  options(scipen=999)
  gtf_HsDf <- as.data.frame(gtfHs)
  RetainedIntron <- gtf_HsDf[gtf_HsDf$transcript_biotype == "retained_intron",]
  RetainedIntron_info <- data.frame(RetainedIntron$ID, unlist(RetainedIntron$Parent), RetainedIntron$transcript_biotype) #latest human annotation (June 2025) required the Parent column to be reformatted to be included in a dataframe

  ExonMeta_filtered <- ExonMeta_filtered[rownames(ExonMeta_filtered) %ni% gsub("_Mega", "", RetainedIntron_info$"unlist.RetainedIntron.Parent."),] #
  IntronMeta_filtered <- IntronMeta_filtered[rownames(IntronMeta_filtered) %ni% gsub("_Mega", "", RetainedIntron_info$"unlist.RetainedIntron.Parent."),] #
  print(identical(rownames(IntronMeta_filtered), rownames(ExonMeta_filtered)))
  
  #rm intron that overlap exon
  gtf_df2 <- as.data.frame(gtf_HsDf)
  gtf_df2_exon <- gtf_df2[gtf_df2$type == "exon",]
  gtf_df2_intron <- gtf_df2[gtf_df2$type == "intron",]
  
  gtf_df2_exon_neg <- gtf_df2_exon[gtf_df2_exon$strand == "-", c(1:3,10)]
  gtf_df2_exon_pos <- gtf_df2_exon[gtf_df2_exon$strand == "+", c(1:3,10)]
  gtf_df2_intron_neg <- gtf_df2_intron[gtf_df2_intron$strand == "-", c(1:3,10)]
  gtf_df2_intron_pos <- gtf_df2_intron[gtf_df2_intron$strand == "+", c(1:3,10)]
  
  overlapping_introns_neg <- list() #
  for (intron in 1:nrow(gtf_df2_intron_neg)) {
    intron_row <- gtf_df2_intron_neg[intron,]
    chromosome <- intron_row[1,1]
    gtf_df2_exon_neg_chromo <- gtf_df2_exon_neg[gtf_df2_exon_neg$seqnames == chromosome,]
    overlap1 <- gtf_df2_exon_neg_chromo[gtf_df2_exon_neg_chromo$start >= intron_row$start & gtf_df2_exon_neg_chromo$end <= intron_row$end,]
    overlap2 <- gtf_df2_exon_neg_chromo[gtf_df2_exon_neg_chromo$start <= intron_row$start & gtf_df2_exon_neg_chromo$end >= intron_row$end,]
    overlap3 <- gtf_df2_exon_neg_chromo[gtf_df2_exon_neg_chromo$start >= intron_row$start & gtf_df2_exon_neg_chromo$end >= intron_row$end & gtf_df2_exon_neg_chromo$start <= intron_row$end,]
    overlap4 <- gtf_df2_exon_neg_chromo[gtf_df2_exon_neg_chromo$start <= intron_row$start & gtf_df2_exon_neg_chromo$end <= intron_row$end & gtf_df2_exon_neg_chromo$end >= intron_row$start,] #& gtf_df2_exon_neg_chromo$end > intron_row$start,]
    if (nrow(overlap1) != 0 | nrow(overlap2) != 0 | nrow(overlap3) != 0 | nrow(overlap4) != 0) {
      overlapping_introns_neg <- append(overlapping_introns_neg, intron_row$ID) #save id ran if overlaps with anything in any of the 4 ways this is possible
    }
  }
  overlapping_introns_pos <- list() #
  for (intron in 1:nrow(gtf_df2_intron_pos)) {
    intron_row <- gtf_df2_intron_pos[intron,]
    chromosome <- intron_row[1,1]
    gtf_df2_exon_pos_chromo <- gtf_df2_exon_pos[gtf_df2_exon_pos$seqnames == chromosome,]
    overlap1 <- gtf_df2_exon_pos_chromo[gtf_df2_exon_pos_chromo$start >= intron_row$start & gtf_df2_exon_pos_chromo$end <= intron_row$end,]
    overlap2 <- gtf_df2_exon_pos_chromo[gtf_df2_exon_pos_chromo$start <= intron_row$start & gtf_df2_exon_pos_chromo$end >= intron_row$end,]
    overlap3 <- gtf_df2_exon_pos_chromo[gtf_df2_exon_pos_chromo$start >= intron_row$start & gtf_df2_exon_pos_chromo$end >= intron_row$end & gtf_df2_exon_pos_chromo$start <= intron_row$end,]
    overlap4 <- gtf_df2_exon_pos_chromo[gtf_df2_exon_pos_chromo$start <= intron_row$start & gtf_df2_exon_pos_chromo$end <= intron_row$end & gtf_df2_exon_pos_chromo$end >= intron_row$start,] #& gtf_df2_exon_pos_chromo$end > intron_row$start,]
    if (nrow(overlap1) != 0 | nrow(overlap2) != 0 | nrow(overlap3) != 0 | nrow(overlap4) != 0) {
      overlapping_introns_pos <- append(overlapping_introns_pos, intron_row$ID)
    }
  }
  
  #get parent ids
  Parents_1 <- gsub("_Mega", "", gtf_HsDf[gtf_HsDf$ID %in% overlapping_introns_neg,]$gene_id)
  Parents_2 <- gsub("_Mega", "", gtf_HsDf[gtf_HsDf$ID %in% overlapping_introns_pos,]$gene_id)
  Parents_all <- c(Parents_1, Parents_2) #

  #remove from intron and exonic
  ExonMeta_filtered_new <- ExonMeta_filtered[rownames(ExonMeta_filtered) %ni% Parents_all,] #   
  IntronMeta_filtered_new <- IntronMeta_filtered[rownames(IntronMeta_filtered) %ni% Parents_all,] #   
  
  #Intron is unspliced
  Unspliced <- IntronMeta_filtered_new
  
  #Exon is spliced
  Spliced <- ExonMeta_filtered_new
  
  print("Introns and exons are in the same order?")
  print(identical(rownames(Spliced), rownames(Unspliced)))

  #normalise and save files
  rownames(Spliced) <- paste0(rownames(Spliced), "-E")
  rownames(Unspliced) <- paste0(rownames(Unspliced), "-I")
  
  png(paste0("RawSpliced_Boxplot_", OutName, ".png"), width=1060, height=680)
  boxplot(log2(Spliced+1), xlab="Sample", ylab="log2(Spliced + 1)")
  dev.off()
  
  png(paste0("RawUnspliced_Boxplot_", OutName, ".png"), width=1060, height=680)
  boxplot(log2(Unspliced+1), xlab="Sample", ylab="log2(Unspliced + 1)")
  dev.off()
  
  if (normalise == "together") {
    Counts <- rbind(Spliced, Unspliced)
    
    png(paste0("RawCounts_Boxplot_", OutName, ".png"), width=1060, height=680)
    boxplot(log2(Counts+1), xlab="Sample", ylab="log2(Counts + 1)")
    dev.off()

    edger <- DGEList(counts=Counts) #first set as DGEList object
    #keep <- rowSums(cpm(edger)>= 5) >= 3
    #edger_filtered <- edger[keep, keep.lib.sizes=FALSE]
    #edger <- calcNormFactors(edger_filtered)
    edger <- calcNormFactors(edger)
    Counts_norm <- log2(cpm(edger, normalized.lib.sizes=TRUE, log=FALSE)+1) #
    Exon_counts_norm <- Counts_norm[grep("-E", rownames(Counts_norm)),] #
    Intron_counts_norm <- Counts_norm[grep("-I", rownames(Counts_norm)),] #
    
    print(identical(nrow(Exon_counts_norm), nrow(Spliced)))
    print(identical(nrow(Intron_counts_norm), nrow(Unspliced)))
    print(identical(nrow(Counts_norm), nrow(Counts)))
    
    print(identical(rownames(Exon_counts_norm), rownames(Spliced)))
    print(identical(rownames(Intron_counts_norm), rownames(Unspliced)))
    
    png(paste0("NormedTogetherCounts_Boxplot_", OutName, ".png"), width=1060, height=680)
    boxplot(Counts_norm, xlab="Sample", ylab="log2(Counts Normalised +1)")
    dev.off() 
    
    png(paste0("NormalisedTogetherUnspliced_Boxplot_", OutName,  ".png"), width=1060, height=680)
    boxplot(Intron_counts_norm, xlab="Sample", ylab="log2(Unspliced Counts Normalised +1)")
    dev.off()
    
    png(paste0("NormalisedTogetherSpliced_Boxplot_", OutName, ".png"), width=1060, height=680)
    boxplot(Exon_counts_norm, xlab="Sample", ylab="log2(Spliced Counts Normalised +1)")
    dev.off()
    
    rownames(Exon_counts_norm) <- gsub("-E", "", rownames(Exon_counts_norm))
    rownames(Intron_counts_norm) <- gsub("-I", "", rownames(Intron_counts_norm))
    
    Exon_counts_norm <- data.frame(rownames(Exon_counts_norm), Exon_counts_norm)
    colnames(Exon_counts_norm)[1] <- "ID"
    
    Intron_counts_norm <- data.frame(rownames(Intron_counts_norm), Intron_counts_norm)
    colnames(Intron_counts_norm)[1] <- "ID"
    
    write.csv(Intron_counts_norm, paste0("Unspliced_NormedTogether_", OutName, ".csv"), row.names=FALSE)
    write.csv(Exon_counts_norm, paste0("Spliced_NormedTogether_", OutName, ".csv"), row.names=FALSE)
  }
  if (normalise == "separately") {
    #norm separately
    edger_exonic <- DGEList(counts=Spliced) #first set as DGEList object
    edger_exonic <- calcNormFactors(edger_exonic)
    Counts_norm_exonic <- log2(cpm(edger_exonic, normalized.lib.sizes=TRUE, log=FALSE)+1) #
    print(identical(nrow(Counts_norm_exonic), nrow(Spliced)))
    
    edger_intronic <- DGEList(counts=Unspliced) #first set as DGEList object
    edger_intronic <- calcNormFactors(edger_intronic)
    Counts_norm_intronic <- log2(cpm(edger_intronic, normalized.lib.sizes=TRUE, log=FALSE)+1) #
    print(identical(nrow(Counts_norm_intronic), nrow(Unspliced)))
    
    png(paste0("NormalisedSeparatelyUnspliced_Boxplot_", OutName, ".png"), width=1060, height=680)
    boxplot(Counts_norm_intronic, xlab="Sample", ylab="log2(Unspliced Counts Normalised +1)")
    dev.off()
    
    png(paste0("NormalisedSeparatelySpliced_Boxplot_", OutName, ".png"), width=1060, height=680)
    boxplot(Counts_norm_exonic, xlab="Sample", ylab="log2(Spliced Counts Normalised +1)")
    dev.off() 
    
    Counts_norm_exonic <- data.frame(rownames(Counts_norm_exonic), Counts_norm_exonic)
    colnames(Counts_norm_exonic)[1] <- "ID"
    
    Counts_norm_intronic <- data.frame(rownames(Counts_norm_intronic), Counts_norm_intronic)
    colnames(Counts_norm_intronic)[1] <- "ID"
    
    Counts_norm_exonic$ID <- gsub("-E", "", Counts_norm_exonic$ID)
    Counts_norm_intronic$ID <- gsub("-I", "", Counts_norm_intronic$ID)
    
    print(identical(Counts_norm_exonic$ID, Counts_norm_intronic$ID))
    print(identical(colnames(Counts_norm_exonic), colnames(Counts_norm_intronic)))
    
    write.csv(Counts_norm_intronic, paste0("Unspliced_NormedSep_", OutName, ".csv"), row.names=FALSE)
    write.csv(Counts_norm_exonic, paste0("Spliced_NormedSep_", OutName, ".csv"), row.names=FALSE)
  }
  print(paste0("Exon-overlapping reads are used as an approximation for spliced transcripts. You have spliced transcript counts for ", nrow(Spliced), " genes and ", ncol(Spliced), " samples."))
  print(paste0("Intron-overlapping reads are used as an approximation for unspliced transcripts. You have unspliced transcript counts for ", nrow(Unspliced), " genes and ", ncol(Unspliced), " samples."))
}

TranscriptInfo <- read.csv("Gene_Transcript_Information.csv") #see the other R script in this sub-directory for how to generate this file from the annotation, this is also given in /Annotation
GeneOccurances <- table(TranscriptInfo$GeneID)
GeneOccurances <- data.frame(names(GeneOccurances), GeneOccurances)
MultipleTrancripts <- GeneOccurances[GeneOccurances$Freq > 1,]
length(unique(TranscriptInfo$TranscriptID))
length(unique(TranscriptInfo$GeneID)) 

MyTranscriptList <- TranscriptInfo[TranscriptInfo$GeneID %in% MultipleTrancripts$Var1,]$TranscriptID

CalculateSplicedUnspliced(Intron="FeatureCounts_IntronMeta.csv", #across-sample compiled gene expression file for the intronic counts
               Exon="FeatureCounts_ExonMeta.csv", #across-sample compiled gene expression file for the exonic counts
               Out="", #specify output directory
               TranscriptList = MyTranscriptList, #see above
               GeneInfo="Gene_Transcript_Information.csv", #see above
               gff3file='Homo_sapiens.GRCh38.104_IntronsFinalSorted.gff3', 'annotation file
               normalise="separately", #normalize introns and exons: "separately" or "together"
               OutName="AllSamples_NormedSeparately") #suffix to add on the end of your spliced and unspliced transcript expression files

CalculateSplicedUnspliced(Intron="FeatureCounts_IntronMeta.csv", #across-sample compiled gene expression file for the intronic counts
              Exon="FeatureCounts_ExonMeta.csv", #across-sample compiled gene expression file for the exonic counts
              Out="", #specify output directory
              TranscriptList = MyTranscriptList, #see above
              GeneInfo="Gene_Transcript_Information.csv", #see above
              gff3file='Homo_sapiens.GRCh38.104_IntronsFinalSorted.gff3', #annotation file
              normalise="together", #normalize introns and exons: "separately" or "together"
              OutName="AllSamples_NormedTogether") #suffix to add on the end of your spliced and unspliced transcript expression files


