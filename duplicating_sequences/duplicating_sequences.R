duplicating_sequences <- function(working_dir,fasta_file_name,freq_file_name) {

# Throwing out error messages if any of the inputs are missing from the command line
x <- 0
error_one <- 0
error_two <- 0
error_three <- 0

killswitch <- "no"

if(missing(working_dir)) {
x <- 1
error_one <- 1
killswitch <- "yes"
}

if(missing(fasta_file_name)) {
x <- 1
error_two <- 2
killswitch <- "yes"
}

if(missing(freq_file_name)) {
x <- 1
error_three <- 3
killswitch <- "yes"
}

if(x==1) {
cat("Call the program by duplicating_sequences(working_dir,fasta_file_name,freq_file_name), where:\nworking_dir == pathway to the folder with your fasta and frequency files e.g. \"C:/blahblahblah\" \nfasta_file_name == the name of your fasta file e.g. \"data.fasta\"\nfreq_file_name == the name of your tab-delimited frequency file, with the sample names \n(including the > fasta prefix) in the first column, and frequency in the second.\n\nExample of input:\nduplicating_sequences(\"C:/Users/Folder/\",\"All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta\",\"frequency.txt\")\n\nSpecific errors/missing inputs:\n")
}
if(error_one==1) {
cat("Sorry, I am missing a working directory pathway\nworking_dir == pathway to the folder with your fasta and frequency files e.g. \"C:/blahblahblah\" \n\n")
}
if(error_two==2) {
cat("Sorry, I am missing a filename for your input fasta file\nfasta_file_name == the name of your fasta file e.g. \"data.fasta\"\n\n")
}

if(error_three==3) {
cat("Sorry, I am missing a filename for your input frequency file\nfreq_file_name == the name of your frequency file e.g. \"frequency.txt\"\n\n")
}


if (killswitch=="yes") {
stop("Please fill in the missing info in your function call")
}

#Checking status of working directory
print(noquote("STEP ONE: Loading in all the variables"))
print(noquote(""))
print(noquote("An error message after this indicates your working directory is not valid"))
flush.console()
setwd(working_dir)
print(noquote("Not to worry, your working directory IS valid! I've successfully set the working directory"))
print(noquote(""))
flush.console()

#Checking status of fasta file
print(noquote("An error message after this indicates your fasta file is not located in the directory you listed"))
flush.console()
intable <- read.table(fasta_file_name,header=FALSE,stringsAsFactors=FALSE,sep="\t")
print(noquote("Not to worry, your fasta file IS located in the directory!"))
print(noquote(""))
flush.console()

#Checking status of fasta file
print(noquote("An error message after this indicates your frequency file is not located in the directory you listed"))
flush.console()
freq <- read.table(freq_file_name,header=FALSE,stringsAsFactors=FALSE,sep="\t")
print(noquote("Not to worry, your frequency file IS located in the directory!"))
print(noquote(""))
flush.console()

rm(error_one)
rm(error_two)
rm(error_three)
rm(killswitch)
rm(fasta_file_name)
rm(freq_file_name)
rm(working_dir)

rows <- dim(intable)[1]

tablelength <- dim(freq)[1]*2

to_write <- matrix(NA,ncol=1,nrow=tablelength)
to_write[1,1] <- paste(freq[1,1],sep="")

to_write_title <- 2
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- paste(freq[(1+(to_write_title/2)),1],sep="")
to_write_title <- to_write_title + 1
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[tablelength,1] <- sequencepaste

rm(rows)
rm(j)
rm(sequencepaste)
rm(to_write_title)
rm(to_write_seq)
rm(x)

intable <- to_write
rm(to_write)

outputlength <- sum(as.numeric(freq[,2]))
output <- matrix(NA, ncol=1,nrow=outputlength*2)
rowno <- 1

for (j in 1:tablelength) {
if (grepl(">",intable[j,1])) {
temp <- intable[j:(j+1),1]

for (k in 1:(tablelength/2)) {
if(intable[j,1]==freq[k,1]) {
upper <- rowno + (freq[k,2]*2) - 1
output[rowno:upper,1] <- temp
rowno <- upper + 1
break
}
}
}
}

write.table(output, "output.fasta",quote=FALSE, col.names=FALSE,row.names=FALSE)
print(noquote("The output of duplicating_sequences has been written to:"))
print("output.fasta")
print(noquote(""))
flush.console()
}
