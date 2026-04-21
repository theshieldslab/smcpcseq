setwd("~/UF/Dental/rob/ana/")

data.ua159 <- read.delim("./UA159/UA159_cds_annotations.tsv", header=F)
colnames(data.ua159) <- c("geneID", "name", "start","end","length","strand","annotation")

origin <- read.delim("./UA159/UA159_TA_Abase_1bp.bed", header=F)
colnames(origin) <- c("insertion","start","end")

### Loop through all insertions
for (i in 1:nrow(origin)) {
  origin[i,1] <- paste0("Insertion_",i)
}
data.ta <- origin

data.ta <- data.frame(start=data.ta[,2], data.ta, geneID="",name="",start=0,end=0,length=0,strand="",annotation="" )

#i.ta=1
#for ( i.ta in 1:nrow(data.ta)) {
#  
#  for ( i.ua159 in 1:nrow(data.ua159)) {
#  
#  if(data.ta[i.ta,1]>=data.ua159[i.ua159,3] & data.ta[i.ta,1]<=data.ua159[i.ua159,4] ){
#    
#    data.ta[i.ta,5] <- data.ua159[i.ua159,1]
#    data.ta[i.ta,6] <- data.ua159[i.ua159,2]
#    data.ta[i.ta,7] <- data.ua159[i.ua159,3]
#    data.ta[i.ta,8] <- data.ua159[i.ua159,4]
#    data.ta[i.ta,9] <- data.ua159[i.ua159,5]
#    data.ta[i.ta,10] <- data.ua159[i.ua159,6]
#    data.ta[i.ta,11] <- data.ua159[i.ua159,7]
#    
#  }
#  
#    
#  }
#}

#write.table(data.ta, "data.ta_annotated.tsv", sep="\t", row.names=F, col.names=T)

### This file is created by lines 17-38. Once created once there is no need to keep running those lines.
data.ta <- read.delim("data.ta_annotated.tsv", header=T)

samples <- read.delim("samples2.txt", header=F)
colnames(samples) <- c("sampleID", "group")

### Merger for counts from all 40 samples
for ( i in 1:40) { 
  
  i.sample <- samples[i,1]
  data.sample <- read.delim(paste0("ta_insertions/",i.sample,"_raw_aln_sorted.TA.depth.tsv"), header=F)
  data.sample <- data.sample[,-1]
  data.sample[,1] <- data.sample[,1]-1
  colnames(data.sample) <- c("start", i.sample)
  
  data.ta <- merge(data.ta, data.sample, by="start")

}
### Removing column 1
data.ta <- data.ta[,-1]

### Writing base dataset 
### This data set has all TA positions gene space + intergenic
### In Excel a new sheet was created containing only gene space insertions
write.xlsx(data.ta, "ta_insertions_annotated_all_samples.xlsx")



