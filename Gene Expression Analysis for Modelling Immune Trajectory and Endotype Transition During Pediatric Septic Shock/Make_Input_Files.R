setwd("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Figure5/")

Spliced <- read.csv("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/Spliced_Unspliced/Spliced_NormedSep_NormSep.csv")
Unspliced <- read.csv("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Counts/compiled/Spliced_Unspliced/Spliced_NormedSep_NormSep.csv")

Endotype <- read.csv("predicted_endotypes_BASICQ_rnaseq.csv")

#Match file and match IDs
LinkerFile <- read.csv("Files_link.csv")
LinkerFile$files <- gsub("~/R Directory/Endotype/BASIC_RNASEQ/basic_quant/", "", LinkerFile$files)
LinkerFile$files <- gsub("_quant/quant.sf", "", LinkerFile$files)

Endotype <- Endotype[Endotype$Sample %in% LinkerFile$X,]
LinkerFile <- LinkerFile[LinkerFile$X %in% Endotype$Sample,]
Endotype <- Endotype[match(LinkerFile$X, Endotype$Sample),]
identical(LinkerFile$X, Endotype$Sample)
Endotype$FileID <- LinkerFile$files
Endotype <- Endotype[Endotype$FileID %in% colnames(Spliced)[2:ncol(Spliced)],]

Person <- gsub(".TP2", "", Endotype$Sample)
Person <- gsub(".TP1", "", Person)
Person <- unique(Person)

Endotype$EndoForFile <- Endotype$Endotype
for (persony in Person) {
  PersonEndoRow <- Endotype[grepl(persony, Endotype$Sample),]
  PersonEndoRow <- PersonEndoRow[order(PersonEndoRow$TP, decreasing=FALSE),]
  PersonEndo <- PersonEndoRow$Endotype
  if (PersonEndo[1] != PersonEndo[2]) {
    Endotype[grepl(persony, Endotype$Sample), 8] <- "E1to2"
  }
}

#Module genes
Blue_Mod <- read.delim("blue_mod_genes.txt")
magenta_Mod <- read.delim("magenta_mod_genes.txt")
Brown_Mod <- read.delim("brown_mod_genes.txt")
red_Mod <- read.delim("red_mod_genes.txt")
turq_Mod <- read.delim("turquoise_mod_genes.txt")

nrow(Blue_Mod) #2805
nrow(magenta_Mod) #143
nrow(Brown_Mod) #1943
nrow(red_Mod) #479
nrow(turq_Mod) #4076

#select genes and generate files

setwd("/media/claired/6be3aab7-4c30-4ac7-9af5-e7a3cc03f466/BASIC/Figure5/Input_Files/")

Sp_blue <- Spliced[Spliced$ID %in% Blue_Mod$gene,]
nrow(Sp_blue) #1579
Unsp_blue <- Unspliced[Unspliced$ID %in% Blue_Mod$gene,]
metadata_blue <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadata_blue) <- c("Sample", "Group")
identical(metadata_blue$Sample, colnames(Sp_blue)[2:ncol(Sp_blue)])
identical(metadata_blue$Sample, colnames(Unsp_blue)[2:ncol(Unsp_blue)])
identical(rownames(Sp_blue), rownames(Unsp_blue))
metadataalt_blue <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadataalt_blue) <- c("Sample", "Group")
identical(metadataalt_blue$Sample, colnames(Sp_blue)[2:ncol(Sp_blue)])
write.csv(Sp_blue, "Spliced_blue.csv", row.names=FALSE)
write.csv(Unsp_blue, "Unspliced_blue.csv", row.names=FALSE)
write.csv(metadata_blue, "Metadata_blue.csv", row.names=FALSE)
write.csv(metadataalt_blue, "MetadataAlt_blue.csv", row.names=FALSE)

Sp_magenta <- Spliced[Spliced$ID %in% magenta_Mod$gene,]
nrow(Sp_magenta) #86
Unsp_magenta <- Unspliced[Unspliced$ID %in% magenta_Mod$gene,]
metadata_magenta <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadata_magenta) <- c("Sample", "Group")
identical(metadata_magenta$Sample, colnames(Sp_magenta)[2:ncol(Sp_magenta)])
identical(metadata_magenta$Sample, colnames(Unsp_magenta)[2:ncol(Unsp_magenta)])
identical(rownames(Sp_magenta), rownames(Unsp_magenta))
metadataalt_magenta <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadataalt_magenta) <- c("Sample", "Group")
identical(metadataalt_magenta$Sample, colnames(Sp_magenta)[2:ncol(Sp_magenta)])
write.csv(Sp_magenta, "Spliced_magenta.csv", row.names=FALSE)
write.csv(Unsp_magenta, "Unspliced_magenta.csv", row.names=FALSE)
write.csv(metadata_magenta, "Metadata_magenta.csv", row.names=FALSE)
write.csv(metadataalt_magenta, "MetadataAlt_magenta.csv", row.names=FALSE)

Sp_brown <- Spliced[Spliced$ID %in% Brown_Mod$gene,]
nrow(Sp_brown) #1079
Unsp_brown <- Unspliced[Unspliced$ID %in% Brown_Mod$gene,]
metadata_brown <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadata_brown) <- c("Sample", "Group")
identical(metadata_brown$Sample, colnames(Sp_brown)[2:ncol(Sp_brown)])
identical(metadata_brown$Sample, colnames(Unsp_brown)[2:ncol(Unsp_brown)])
identical(rownames(Sp_brown), rownames(Unsp_brown))
metadataalt_brown <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadataalt_brown) <- c("Sample", "Group")
identical(metadataalt_brown$Sample, colnames(Sp_brown)[2:ncol(Sp_brown)])
write.csv(Sp_brown, "Spliced_brown.csv", row.names=FALSE)
write.csv(Unsp_brown, "Unspliced_brown.csv", row.names=FALSE)
write.csv(metadata_brown, "Metadata_brown.csv", row.names=FALSE)
write.csv(metadataalt_brown, "MetadataAlt_brown.csv", row.names=FALSE)

Sp_red <- Spliced[Spliced$ID %in% red_Mod$gene,]
nrow(Sp_red) #300
Unsp_red <- Unspliced[Unspliced$ID %in% red_Mod$gene,]
metadata_red <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadata_red) <- c("Sample", "Group")
identical(metadata_red$Sample, colnames(Sp_red)[2:ncol(Sp_red)])
identical(metadata_red$Sample, colnames(Unsp_red)[2:ncol(Unsp_red)])
identical(rownames(Sp_red), rownames(Unsp_red))
  metadataalt_red <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadataalt_red) <- c("Sample", "Group")
identical(metadataalt_red$Sample, colnames(Sp_red)[2:ncol(Sp_red)])
write.csv(Sp_red, "Spliced_red.csv", row.names=FALSE)
write.csv(Unsp_red, "Unspliced_red.csv", row.names=FALSE)
write.csv(metadata_red, "Metadata_red.csv", row.names=FALSE)
write.csv(metadataalt_red, "MetadataAlt_red.csv", row.names=FALSE)

Sp_blue <- Spliced[Spliced$ID %in% Blue_Mod$gene,]
nrow(Sp_blue) #1579
Unsp_blue <- Unspliced[Unspliced$ID %in% Blue_Mod$gene,]
metadata_blue <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadata_blue) <- c("Sample", "Group")
identical(metadata_blue$Sample, colnames(Sp_blue)[2:ncol(Sp_blue)])
identical(metadata_blue$Sample, colnames(Unsp_blue)[2:ncol(Unsp_blue)])
identical(rownames(Sp_blue), rownames(Unsp_blue))
metadataalt_blue <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadataalt_blue) <- c("Sample", "Group")
identical(metadataalt_blue$Sample, colnames(Sp_blue)[2:ncol(Sp_blue)])
write.csv(Sp_blue, "Spliced_blueTPCol.csv", row.names=FALSE)
write.csv(Unsp_blue, "Unspliced_blueTPCol.csv", row.names=FALSE)
write.csv(metadata_blue, "Metadata_blueTPCol.csv", row.names=FALSE)
write.csv(metadataalt_blue, "MetadataAlt_blueTPCol.csv", row.names=FALSE)

Sp_magenta <- Spliced[Spliced$ID %in% magenta_Mod$gene,]
nrow(Sp_magenta) #86
Unsp_magenta <- Unspliced[Unspliced$ID %in% magenta_Mod$gene,]
metadata_magenta <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadata_magenta) <- c("Sample", "Group")
identical(metadata_magenta$Sample, colnames(Sp_magenta)[2:ncol(Sp_magenta)])
identical(metadata_magenta$Sample, colnames(Unsp_magenta)[2:ncol(Unsp_magenta)])
identical(rownames(Sp_magenta), rownames(Unsp_magenta))
metadataalt_magenta <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadataalt_magenta) <- c("Sample", "Group")
identical(metadataalt_magenta$Sample, colnames(Sp_magenta)[2:ncol(Sp_magenta)])
write.csv(Sp_magenta, "Spliced_magentaTPCol.csv", row.names=FALSE)
write.csv(Unsp_magenta, "Unspliced_magentaTPCol.csv", row.names=FALSE)
write.csv(metadata_magenta, "Metadata_magentaTPCol.csv", row.names=FALSE)
write.csv(metadataalt_magenta, "MetadataAlt_magentaTPCol.csv", row.names=FALSE)

Sp_brown <- Spliced[Spliced$ID %in% Brown_Mod$gene,]
nrow(Sp_brown) #1079
Unsp_brown <- Unspliced[Unspliced$ID %in% Brown_Mod$gene,]
metadata_brown <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadata_brown) <- c("Sample", "Group")
identical(metadata_brown$Sample, colnames(Sp_brown)[2:ncol(Sp_brown)])
identical(metadata_brown$Sample, colnames(Unsp_brown)[2:ncol(Unsp_brown)])
identical(rownames(Sp_brown), rownames(Unsp_brown))
metadataalt_brown <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadataalt_brown) <- c("Sample", "Group")
identical(metadataalt_brown$Sample, colnames(Sp_brown)[2:ncol(Sp_brown)])
write.csv(Sp_brown, "Spliced_brownTPCol.csv", row.names=FALSE)
write.csv(Unsp_brown, "Unspliced_brownTPCol.csv", row.names=FALSE)
write.csv(metadata_brown, "Metadata_brownTPCol.csv", row.names=FALSE)
write.csv(metadataalt_brown, "MetadataAlt_brownTPCol.csv", row.names=FALSE)

Sp_red <- Spliced[Spliced$ID %in% red_Mod$gene,]
nrow(Sp_red) #300
Unsp_red <- Unspliced[Unspliced$ID %in% red_Mod$gene,]
metadata_red <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadata_red) <- c("Sample", "Group")
identical(metadata_red$Sample, colnames(Sp_red)[2:ncol(Sp_red)])
identical(metadata_red$Sample, colnames(Unsp_red)[2:ncol(Unsp_red)])
identical(rownames(Sp_red), rownames(Unsp_red))
metadataalt_red <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadataalt_red) <- c("Sample", "Group")
identical(metadataalt_red$Sample, colnames(Sp_red)[2:ncol(Sp_red)])
write.csv(Sp_red, "Spliced_redTPCol.csv", row.names=FALSE)
write.csv(Unsp_red, "Unspliced_redTPCol.csv", row.names=FALSE)
write.csv(metadata_red, "Metadata_redTPCol.csv", row.names=FALSE)
write.csv(metadataalt_red, "MetadataAlt_redTPCol.csv", row.names=FALSE)

Sp_turq <- Spliced[Spliced$ID %in% turq_Mod$gene,]
nrow(Sp_turq) #2589
Unsp_turq <- Unspliced[Unspliced$ID %in% turq_Mod$gene,]
metadata_turq <- data.frame(Endotype$FileID, Endotype$EndoForFile)
colnames(metadata_turq) <- c("Sample", "Group")
identical(metadata_turq$Sample, colnames(Sp_turq)[2:ncol(Sp_turq)])
identical(metadata_turq$Sample, colnames(Unsp_turq)[2:ncol(Unsp_turq)])
identical(rownames(Sp_turq), rownames(Unsp_turq))
metadataalt_turq <- data.frame(Endotype$FileID, Endotype$TP)
colnames(metadataalt_turq) <- c("Sample", "Group")
identical(metadataalt_turq$Sample, colnames(Sp_turq)[2:ncol(Sp_turq)])
write.csv(Sp_turq, "Spliced_turq.csv", row.names=FALSE)
write.csv(Unsp_turq, "Unspliced_turq.csv", row.names=FALSE)
write.csv(metadata_turq, "Metadata_turq.csv", row.names=FALSE)
write.csv(metadataalt_turq, "MetadataAlt_turq.csv", row.names=FALSE)


blue <- read.csv("Spliced_magenta.csv")
magenta <- read.csv("Spliced_blue.csv")
brown <- read.csv("Spliced_brown.csv")
red <- read.csv("Spliced_red.csv")
turquoise <- read.csv("Spliced_turq.csv")

blue[blue$ID == "ENSG00000091317",] #no
magenta[magenta$ID == "ENSG00000146425",] #no
red[red$ID == "ENSG00000129472",] #yes
brown[brown$ID == "ENSG00000171791",] #no
turquoise[turquoise$ID == "ENSG00000132341",] #YES




