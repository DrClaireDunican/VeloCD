Exon_Count_files <- list.files(path = "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/", pattern = "_ExonMeta.txt")

`%ni%` = Negate(`%in%`)
summary_files <- grep('summary', Exon_Count_files, value=TRUE)
Exon_Count_files <- Exon_Count_files[Exon_Count_files %ni% summary_files]

length(Exon_Count_files) #144

setwd("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/")

#Exons:
ExonFile <- read.delim('IGF126586_S1_L001_R1_001.fastq.gz_ExonMeta.txt', comment.char = '#') 
ExonFileUseful <- ExonFile[,c(1,6:7)]
colnames(ExonFileUseful) <- c(colnames(ExonFileUseful)[1:2], "Count")
Exon_count_DataFrame <- ExonFileUseful
write.csv(ExonFileUseful[,1:2], "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/Lengths_ExonMeta.csv", row.names=FALSE)
colnames(Exon_count_DataFrame) <- c("ID", "Length", "IGF126586_S1_L001_R1_001")
Exon_count_DataFrame <- Exon_count_DataFrame[,-2]
ExonFileUseful <- NULL
for (Exon_Count_file in Exon_Count_files[2:length(Exon_Count_files)]) {
  ExonFile <- read.delim(Exon_Count_file, comment.char = '#')
  ExonFileUseful <- ExonFile[,c(1,7)]
  print(Exon_Count_file)
  print(dim(ExonFileUseful))
  print(identical(ExonFileUseful$Geneid, Exon_count_DataFrame$ID)) #check the Exons are in the same order
  colnames(ExonFileUseful) <- c("ID", "Count")
  Exon_count_DataFrame <- data.frame(Exon_count_DataFrame, ExonFileUseful$Count)
  new_col_name <- toString(Exon_Count_file)
  new_col_name <- gsub('.fastq.gz_ExonMeta.txt', "", new_col_name)
  print(new_col_name)
  colnames(Exon_count_DataFrame) <- c(colnames(Exon_count_DataFrame)[1:ncol(Exon_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
} #60664 2

write.csv(Exon_count_DataFrame, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_ExonMeta.csv", row.names=FALSE)
dim(Exon_count_DataFrame) #60672 81
