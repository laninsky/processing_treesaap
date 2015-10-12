#Example usage
#haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil", "output.fasta")
#haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil/test", "All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta")


#working_dir <- "C:/Users/Alana/Dropbox/Mitogenome_Phil"
#file_name <- "output.fasta"
#file_name <- "All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta"


haplotyping_sequences <- function(working_dir,file_name) {
# The stringr library is required
library(stringr)

# Throwing out error messages if any of the inputs are missing from the command line
x <- 0
error_one <- 0
error_two <- 0

killswitch <- "no"

if(missing(working_dir)) {
x <- 1
error_one <- 1
killswitch <- "yes"
}

if(missing(file_name)) {
x <- 1
error_two <- 2
killswitch <- "yes"
}

if(x==1) {
cat("Call the program by haplotyping_sequences(working_dir,file_name), where:\nworking_dir == pathway to the folder with your fasta file e.g. \"C:/blahblahblah\" \nfile_name == the name of your fasta file e.g. \"data.fas\"\n\nExample of input:\nhaplotyping_sequences(\"C:/Users/Folder/\",\"ATL_by_region_394.fasta\")\n\nSpecific errors/missing inputs:\n")
}
if(error_one==1) {
cat("Sorry, I am missing a working directory pathway\nworking_dir == pathway to the folder with your fasta file e.g. \"C:/blahblahblah\" \n\n")
}
if(error_two==2) {
cat("Sorry, I am missing a filename for your input fasta file\nfile_name == the name of your fasta file e.g. \"data.fasta\"\n\n")
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
intable <- read.table(file_name,header=FALSE,stringsAsFactors=FALSE,sep="\t")
print(noquote("Not to worry, your fasta file IS located in the directory!"))
print(noquote(""))
flush.console()

rm(error_one)
rm(error_two)
rm(killswitch)
rm(file_name)
rm(working_dir)
rm(x)

#If there are line-breaks, getting all the sequence on to a single line, with sequence name in col 1, sequence in col 2

rows <- dim(intable)[1]

haplist <- t(as.matrix(c(intable[1,1],NA)))
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
temp <- t(as.matrix(c(intable[j,1],NA)))
haplist <- rbind(haplist,temp)
haplist_rows <- dim(haplist)[1]
haplist[(haplist_rows-1),2] <- sequencepaste
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

haplist_rows <- dim(haplist)[1]
haplist[haplist_rows,2] <- sequencepaste

rm(rows)
rm(j)
rm(temp)
rm(haplist_rows)
rm(sequencepaste)
rm(intable)

haplistlength <- dim(haplist)[1]

checkingforduphapnames <- as.data.frame(table(haplist[2:haplistlength,1]))
checkingforduphapnames <- checkingforduphapnames[order(checkingforduphapnames[,2]),]
checkingforduphapnameslen <- dim(checkingforduphapnames)[1]
if(checkingforduphapnames[checkingforduphapnameslen,2]>1) {
print(noquote(""))
print(noquote("Multiple haplotype definitions with the same haplotype name are present"))
print(noquote("Please check the following haplotype name, and make sure all haplotype definitions"))
print(noquote("have unique haplotype names:"))
print(droplevels(checkingforduphapnames[checkingforduphapnameslen,1]))
flush.console()
stop("Please rename the multiple haplotype definitions that have the above name so that they have unique names")
}

minhaplength <- min(nchar(haplist[2:haplistlength,2]))
maxhaplength <- max(nchar(haplist[2:haplistlength,2]))

#Error messages tested and funcitons
if(minhaplength==maxhaplength) {
   if(all(is.na(haplist[2:haplistlength,2]))) {
stop("\n\nNo sequence given for haplotypes. Please make sure you have sequence defined for each of your haplotypes\n")
}
} else  {
stop("\n\nYour file does not have equal length haplotypes. Please check your alignment and try again\n")
}

print(noquote("Your haplotypes appear to be of equal lengths"))
print(noquote(""))
flush.console()

# Getting some parameters that will be used whether haplotype and/or nucleotide diversity is being tested

no_seqs <- dim(haplist)[1]  #Getting the number of haplotypes from the input file

print(noquote("I am now going to check for duplicate haplotypes in your haplotype definitions"))
flush.console()

# We need to check if any of the haplotypes are identical and pool these if so
# Defining our variables for the loop below
pattern <- NULL
propdiffs <- matrix(NA,nrow = (no_seqs + 1),ncol = (no_seqs + 1))
propdiffs[1,2:(no_seqs+1)] <- t(haplist[,1])
propdiffs[2:(no_seqs+1),1] <- haplist[,1]
ambigs <- c("R","Y","S","W","K","M","B","D","H","V","N","r","y","s","w","k","m","b","d","h","v","n")
OKbases <- c("A","C","G","T","-","a","c","g","t")

# Calculating the proportion of differences between haplotypes for all base positions that do not have an ambiguous nucleotide
print(noquote("Will update you on progress for calculating proportion of differences between sequences at every fifth sequence"))
print(noquote("Currently up to the following sequence in your file"))
flush.console()

for (j in 1:no_seqs) {

if(round(j/5)==(j/5)) {
print(noquote(j))
flush.console()
}

k <- j + 1

first <- unlist(strsplit(haplist[j,2],pattern))
no_bp <- length(first)
while (k <= no_seqs) {
if(!(haplist[j,2]==haplist[k,2])) {
tot_bp <- 0
mismatch <- 0
second <- unlist(strsplit(haplist[k,2],pattern))
for (m in 1:no_bp) {
if (!(first[m] %in% ambigs)) {
if (!(first[m] %in% OKbases)) {
stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
}
if (!(second[m] %in% ambigs)) {
if (!(second[m] %in% OKbases)) {
stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
}

tot_bp <- tot_bp + 1

if (first[m]!=second[m]) {
if (tolower(first[m])==second[m]) {
break
}
if (toupper(first[m])==second[m]) {
break
}
mismatch <- mismatch + 1 }
                             }
                           }
                    }
propdiffs[(k+1),(j+1)] <- mismatch/tot_bp
                            } else {
propdiffs[(k+1),(j+1)] <- 0
}

k <- k + 1

                          }
}

rm(ambigs)
rm(checkingforduphapnameslen)
rm(checkingforduphapnames)
rm(j)
rm(k)
rm(m)
rm(maxhaplength)
rm(minhaplength)
rm(mismatch)
rm(no_bp)
rm(first)
rm(OKbases)
rm(pattern)
rm(second)
rm(tot_bp)
rm(haplistlength)

j <- 1
namearray <- NULL

while (j <= no_seqs) {
namearray[j] <- haplist[j,1]
k <- j + 1
m <- 0
while (k <= no_seqs) {
if (propdiffs[(k+1),(j+1)]==0) {
m <- m + 1
namearray[j] <- paste(namearray[j],haplist[k,1])
}
k <- k + 1
}
j <- j + 1
}

finnamearray <- NULL

namearray <- cbind(namearray,NA)
namearraylen <- dim(namearray)[1]

j <- 1
i <- 1

if(namearraylen==1) {
finnamearray <- namearray
} else {
for (j in 1:(namearraylen)) {
if(is.na(namearray[j,2])) {
namearray[j,2] <- i
rowlen <- length(unlist(strsplit(namearray[j,1], " ")))
for (m in 1:rowlen) {
k <- j + 1
while (k <= namearraylen) {
if(is.na(namearray[k,2])) {
if(((unlist(strsplit(namearray[j,1]," "))) [m]) %in% (unlist(strsplit(namearray[k,1]," ")))) {
namearray[k,2] <- i
}
}
k <- k + 1
}
}
i <- i + 1
}
j <- j + 1
}

if(is.na(namearray[namearraylen,2])) {
namearray[namearraylen,2] <- i + 1
}

namearray <- namearray[order(namearray[,2]),]
i <- 1
j <- 1

while (j < namearraylen) {
k <- j + 1
finnamearray[i] <- namearray[j,1]
while ((k <= namearraylen) && (namearray[j,2]==namearray[k,2])) {
finnamearray[i] <- paste(finnamearray[i],namearray[k,1])
k <- k + 1
}
i <- i + 1
j <- k
}
}

rm(j)
rm(k)
rm(m)
rm(namearray)
rm(namearraylen)

uniques_for_adding <- NULL
sets_of_unique_haplotypes <- "sets of unique haplotypes, separated on new lines"
j <- 1

while (j <= length(finnamearray)) {
vectorizing <- unlist(strsplit(finnamearray[j], " "))
onlyunique <- vectorizing[!duplicated(vectorizing)]
if(length(onlyunique)>1) {
print(noquote("The following set of haplotypes is identical"))
print(noquote(onlyunique))
print(noquote(""))
flush.console()
}
uniques_for_adding <- rbind(uniques_for_adding,onlyunique[1])

for (m in 1:no_seqs) {
if (haplist[m,1]==onlyunique[1]) {
uniques_for_adding <- rbind(uniques_for_adding,haplist[m,2])
}
}
sets_of_unique_haplotypes <- rbind(sets_of_unique_haplotypes,(paste0(onlyunique,collapse=", ")))
j <- j + 1
}

print(noquote("The fasta file of the unique haplotypes has been written to:"))
print("haplist.fasta")
print(noquote(""))
flush.console()
write.table(uniques_for_adding, "haplist.fasta", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)

print(noquote("The sets of duplicate sequences for each unique haplotypes has been written to:"))
print("duplicate_sequences.txt")
print(noquote(""))
flush.console()
write.table(sets_of_unique_haplotypes, "duplicate_sequences.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)
}
