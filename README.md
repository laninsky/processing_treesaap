# processing_treesaap
Summarizing results over multiple pairwise comparisons using R

```
#Set your working directory in the folder containing all the subfolders for the comparison
library(stringr)
setwd("C:/Users/a499a400/Dropbox/Mitogenome_Phil/selection_analyses/treesaap/multi_pma/pairwise/")

filelist <- list.files()

#creating list of files which have Pma-Pma comparisons and excluding them (Pmas are all prefixed by mtgen)
newfilelist <- filelist[which((str_count(filelist[],"MTGEN"))<2)]
nocomps <- length(filelist)
output <- c("codon","comparison","property")
write.table(output, "output.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

for (i in 1:nocomps) {
toprint <- paste(round((i/nocomps*100),3),"% through the files",sep="")
print(toprint)
flush.console()
tempcodons <- readLines(paste(newfilelist[i],"/Evpthwy/Subs-SigProp.txt",sep=""))
for (j in 2:(length(tempcodons))) {
templine <- unlist(strsplit(tempcodons[j],"\t"))
propertylists <- which(grepl("\\(6, ", templine[]) | grepl("\\(7, ", templine[]) | grepl("\\(8, ", templine[]) )
if(length(propertylists)>0) {
for (k in 1:length(propertylists)) {
tempproperty <- unlist(strsplit(templine[propertylists[k]],"\\) "))[2]
tempoutput <- c(templine[1],templine[2],tempproperty)
write.table(tempoutput, "output.txt",append=TRUE,quote=FALSE, col.names=FALSE,row.names=FALSE)
}
}
}
}


