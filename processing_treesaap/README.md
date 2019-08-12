# processing_treesaap
Summarizing results over multiple pairwise comparisons using R

```
#Set your working directory in the folder containing all the subfolders for the comparison
library(stringr)
library(data.table)
library(plyr)

setwd("C:/Users/a499a400/Dropbox/Mitogenome_Phil/selection_analyses/treesaap/multi_pma/pairwise/")

filelist <- list.files()

#creating list of files which have Pma comparisons and excluding them (Pmas are all prefixed by mtgen)
#newfilelist <- filelist[which((str_count(filelist[],"MTGEN"))<2)]


#creating list of files which have a Pma comparison and including them (Pmas are all prefixed by mtgen)
newfilelist <- filelist[which((str_count(filelist[],"MTGEN"))==1)]

nocomps <- length(newfilelist)
output <- matrix(c("codon","comparison","property"),nrow=1)
write.table(output, "output.txt",quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")


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
tempoutput <- matrix(c(templine[1],templine[2],tempproperty),nrow=1)
write.table(tempoutput, "output.txt",quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
}
}
}
}

tempoutput <- read.table("output.txt",sep="\t",header=TRUE)
outputtable <- as.data.table(tempoutput)
out <- ddply(outputtable, .(codon,property), as.data.frame(nrow))
out <- out[(which(out[,3]==nocomps)),]
```

### Version history
This script was written for analyses associated with:  
Morin, P.A., Foote, A.D., Baker, C.S., Hancock‐Hanser, B.L., Kaschner, K., Mate, B.R., Mesnick, S.L., Pease, V.L., Rosel, P.E. and Alexander, A., 2018. Demography or selection on linked cultural traits or genes? Investigating the driver of low mtDNA diversity in the sperm whale using complementary mitochondrial and nuclear genome analyses. Molecular ecology, 27(11), pp.2604-2619.
