library(conumee)
library(minfi); library(minfiData)
COSMIC <- read.csv("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/COSMIC_Cancer_Genes/Cosmic Cancer Gene Census_all hg37.tsv", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

library(dplyr)

COSMIC <- COSMIC %>%
  select("Gene Symbol", "Genome Location", "Chr Band", "Entrez GeneId",   "Name") %>%
  mutate(chr = sapply(strsplit(COSMIC$`Genome Location`, "\\:|-"), '[[', 1), start = sapply(strsplit(COSMIC$`Genome Location`, "\\:|-"), '[[', 2), end = sapply(strsplit(COSMIC$`Genome Location`, "\\:|-"), '[[', 3)) %>%
  rename(name="Gene Symbol", location="Genome Location", EntrezID= "Entrez GeneId",band = "Chr Band", geneName= "Name") %>%
  #select(name, start, end, chr,EntrezID, band) %>%
  select(name, start, end, chr) %>%
  mutate(chr=paste0("chr", chr)) 
#COSMIC <- COSMIC[COSMIC$name %in%probe.features$gene,]## Out of 717 only 703 have annotation in probe.features of minfi.


library(GenomicRanges)
COSMIC_regions <- makeGRangesFromDataFrame(COSMIC, keep.extra.columns = TRUE)
anno <- CNV.create_anno(detail_regions = COSMIC_regions, chrXY = TRUE)
COSMIC_regions


#COSMIC.1 <- COSMIC %>%
#  select(chr, start, end, name)
#write.csv(COSMIC.1, file="COSMIC_1.bed", row.names=FALSE) ## make it tab delimited file without header line
#OSMIC_regions <- import("COSMIC_1.bed", format = "BED")
#anno <- CNV.create_anno(detail_regions = COSMIC_regions, chrXY = TRUE)


data(MsetEx)
d <- CNV.load(MsetEx)
x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'], 
                                            ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))

CNV.genomeplot(x)
CNV.detailplot(x, name="CYSLTR2")

head(CNV.write(x, what='detail'))
