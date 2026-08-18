Gene_Count_files <- list.files(path = "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/", pattern = "_Gene.txt")

`%ni%` = Negate(`%in%`)
summary_files <- grep('summary', Gene_Count_files, value=TRUE)
Gene_Count_files <- Gene_Count_files[Gene_Count_files %ni% summary_files]

length(Gene_Count_files) #80

setwd("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/")

#Genes:
GeneFile <- read.delim('IGF126586_S1_L001_R1_001.fastq.gz_Gene.txt', comment.char = '#') 
GeneFileUseful <- GeneFile[,c(1,6:7)]
colnames(GeneFileUseful) <- c(colnames(GeneFileUseful)[1:2], "Count")
Gene_count_DataFrame <- GeneFileUseful
write.csv(GeneFileUseful[,1:2], "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/Lengths_Gene.csv", row.names=FALSE)
colnames(Gene_count_DataFrame) <- c("ID", "Length", "IGF126586_S1_L001_R1_001")
Gene_count_DataFrame <- Gene_count_DataFrame[,-2]
GeneFileUseful <- NULL
for (Gene_Count_file in Gene_Count_files[2:length(Gene_Count_files)]) {
  GeneFile <- read.delim(Gene_Count_file, comment.char = '#')
  GeneFileUseful <- GeneFile[,c(1,7)]
  print(Gene_Count_file)
  print(dim(GeneFileUseful))
  print(identical(GeneFileUseful$Geneid, Gene_count_DataFrame$ID)) #check the Genes are in the same order
  colnames(GeneFileUseful) <- c("ID", "Count")
  Gene_count_DataFrame <- data.frame(Gene_count_DataFrame, GeneFileUseful$Count)
  new_col_name <- toString(Gene_Count_file)
  new_col_name <- gsub('.fastq.gz_Gene.txt', "", new_col_name)
  print(new_col_name)
  colnames(Gene_count_DataFrame) <- c(colnames(Gene_count_DataFrame)[1:ncol(Gene_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
} #60664 2

write.csv(Gene_count_DataFrame, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_Gene.csv", row.names=FALSE)
dim(Gene_count_DataFrame) #60664 81

