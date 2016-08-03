# processing_treesaap
Summarizing results over multiple pairwise comparisons using R

```
#Set your working directory in the folder containing all the subfolders for the comparison
setwd("C:/Users/a499a400/Dropbox/Mitogenome_Phil/selection_analyses/treesaap/multi_pma/pairwise/")

filelist <- list.files()

#creating list of files which have Pma-Pma comparisons and excluding them (Pmas are all prefixed by mtgen)
newfilelist <- filelist[which((str_count(filelist[],"MTGEN"))<2)]
nocomps <- length(filelist)

for (i in 1:nocomps) {
tempcodons <- read.table(paste(newfilelist[i],"/Evpthwy/SigCodons.txt",sep=""),sep="\t",skip=1)
tempproperties <- 
