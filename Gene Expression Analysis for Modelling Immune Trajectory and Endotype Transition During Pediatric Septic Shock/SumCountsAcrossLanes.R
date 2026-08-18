Gene_count_DataFrame <- read.csv("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_Gene.csv")
Exon_count_DataFrame <- read.csv("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_ExonMeta.csv")
Intron_count_DataFrame <- read.csv("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_IntronMeta.csv")

SampleNames <- unique(gsub("_L.*", "", colnames(Gene_count_DataFrame)[2:ncol(Gene_count_DataFrame)])) #40

Gene_Summed <- data.frame(Gene_count_DataFrame$ID)
Exon_Summed <- data.frame(Exon_count_DataFrame$ID)
Intron_Summed <- data.frame(Intron_count_DataFrame$ID)

for (sampleName in SampleNames) {
  #Genes
  MySampleCols_Gene <- Gene_count_DataFrame[, grepl(sampleName, colnames(Gene_count_DataFrame))]
  NewSumCol_Gene <- rowSums(MySampleCols_Gene)
  Gene_Summed <- data.frame(Gene_Summed, NewSumCol_Gene)
  colnames(Gene_Summed)[ncol(Gene_Summed)] <- sampleName
  #Exons
  MySampleCols_Exon <- Exon_count_DataFrame[, grepl(sampleName, colnames(Exon_count_DataFrame))]
  NewSumCol_Exon <- rowSums(MySampleCols_Exon)
  Exon_Summed <- data.frame(Exon_Summed, NewSumCol_Exon)
  colnames(Exon_Summed)[ncol(Exon_Summed)] <- sampleName
  #Intron
  MySampleCols_Intron <- Intron_count_DataFrame[, grepl(sampleName, colnames(Intron_count_DataFrame))]
  NewSumCol_Intron <- rowSums(MySampleCols_Intron)
  Intron_Summed <- data.frame(Intron_Summed, NewSumCol_Intron)
  colnames(Intron_Summed)[ncol(Intron_Summed)] <- sampleName
}

colnames(Gene_Summed)[1] <- "ID"
write.csv(Gene_Summed, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_Gene_Summed.csv", row.names=FALSE)

colnames(Exon_Summed)[1] <- "ID"
write.csv(Exon_Summed, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_Exon_Summed.csv", row.names=FALSE)

colnames(Intron_Summed)[1] <- "ID"
write.csv(Intron_Summed, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_Intron_Summed.csv", row.names=FALSE)


