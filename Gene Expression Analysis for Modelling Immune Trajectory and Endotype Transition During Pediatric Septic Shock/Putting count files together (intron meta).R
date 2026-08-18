Intron_Count_files <- list.files(path = "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/", pattern = "_IntronMeta.txt")

`%ni%` = Negate(`%in%`)
summary_files <- grep('summary', Intron_Count_files, value=TRUE)
Intron_Count_files <- Intron_Count_files[Intron_Count_files %ni% summary_files]

length(Intron_Count_files) #80

setwd("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/")

#Intron_s:
Intron_File <- read.delim('IGF126586_S1_L001_R1_001.fastq.gz_IntronMeta.txt', comment.char = '#') 
Intron_FileUseful <- Intron_File[,c(1,6:7)]
colnames(Intron_FileUseful) <- c(colnames(Intron_FileUseful)[1:2], "Count")
Intron_count_DataFrame <- Intron_FileUseful
write.csv(Intron_FileUseful[,1:2], "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/Lengths_IntronMeta.csv", row.names=FALSE)
colnames(Intron_count_DataFrame) <- c("ID", "Length", "IGF126586_S1_L001_R1_001")
Intron_count_DataFrame <- Intron_count_DataFrame[,-2]
Intron_FileUseful <- NULL
for (Intron_Count_file in Intron_Count_files[2:length(Intron_Count_files)]) {
  Intron_File <- read.delim(Intron_Count_file, comment.char = '#')
  Intron_FileUseful <- Intron_File[,c(1,7)]
  print(Intron_Count_file)
  print(dim(Intron_FileUseful))
  print(identical(Intron_FileUseful$Geneid, Intron_count_DataFrame$ID)) #check the Intron_s are in the same order
  colnames(Intron_FileUseful) <- c("ID", "Count")
  Intron_count_DataFrame <- data.frame(Intron_count_DataFrame, Intron_FileUseful$Count)
  new_col_name <- toString(Intron_Count_file)
  new_col_name <- gsub('.fastq.gz_IntronMeta.txt', "", new_col_name)
  print(new_col_name)
  colnames(Intron_count_DataFrame) <- c(colnames(Intron_count_DataFrame)[1:ncol(Intron_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
} #37683 2

write.csv(Intron_count_DataFrame, "/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/FeatureCounts_IntronMeta.csv", row.names=FALSE)
dim(Intron_count_DataFrame) #37683   81
