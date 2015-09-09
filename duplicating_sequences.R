duplicating_sequences <- function(working_dir,fasta_file_name,freq_file_name) {

# The stringr library is required
library(stringr)

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

#Checking status of arlequin file
print(noquote("An error message after this indicates your file is not located in the directory you listed"))
flush.console()
input <- readLines(file_name)
print(noquote("Not to worry, your file IS located in the directory! I'm pulling it into my memory to extract the parts we are interested in"))
print(noquote(""))
flush.console()
