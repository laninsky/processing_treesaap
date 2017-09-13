working_dir <- "C:\\Users\\a499a400\\Dropbox\\Mitogenome_Phil"
file <- "clade_frequencies.txt.R"
no_permuts <- 10

subsampling_clades <- function(working_dir,file,no_permuts) {

#checking to make sure that all arguments are present
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

if(missing(file)) {
x <- 1
error_two <- 2
killswitch <- "yes"
}

if(missing(no_permuts)) {
x <- 1
error_three <- 3
killswitch <- "yes"
}

if(x==1) {
cat("Call the program by subsampling_clades(working_dir,file,no_permuts), where:\nworking_dir == pathway to the folder with your clade file e.g. \"C:/blahblahblah\" \nfile == the name of your clade file (described in readme) e.g. \"data.txt\"\nno_permuts == the number of permutations you wish to perform e.g. 1000\n\nExample of input:\nsubsampling_clades(\"C:/Users/Folder/\",\"ATL_by_region_394.txt\",1000)\n\nSpecific errors/missing inputs:\n")
}
if(error_one==1) {
cat("Sorry, I am missing a working directory pathway\nworking_dir == pathway to the folder with your clade file e.g. \"C:/blahblahblah\" \n\n")
}
if(error_two==2) {
cat("Sorry, I am missing a filename for your clade file\nfile_name == the name of your clade file e.g. \"data.txt\"\n\n")
}
if(error_three==3) {
cat("How many permutations would you like to complete?\nno_permuts == the number of iterations you wish to perform e.g. 1000\n\n")
}

# Error message tested and functioning
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

#Checking file
print(noquote("An error message after this indicates your file is not located in the directory you listed"))
flush.console()
input <- as.matrix(read.table(file))
print(noquote("Not to worry, your file IS located in the directory!"))
print(noquote(""))
flush.console()

#Getting_sample_sizes_for_each_region
samplesizes <- matrix(NA,ncol=2,nrow=(dim(input)[2]-1))
samplesizes[,1] <- input[1,2:(dim(input)[2])]
samplesizes[,2] <- colSums(matrix(as.numeric(input[2:(dim(input)[1]),2:(dim(input)[2])]),ncol=(dim(input)[2]-1)))
samplesizes <- samplesizes[order(as.numeric(samplesizes[,2])),]

for (i in 1:(dim(samplesizes)[1]-1)) {
   record_matrix <- NULL
   for (j in 2:(dim(input)[2])) {
       if ((as.numeric(samplesizes[which(input[1,j]==samplesizes[,1]),2]))>as.numeric(samplesizes[i,2])) {
          temp_array <- NULL
          for (k in 2:(dim(input)[1])) {
            if (as.numeric(input[k,j])>0) {
              temp_array <- c(temp_array,rep(input[k,1],as.numeric(input[k,j])))
            }
          }  
          temp_matrix <- matrix(NA,nrow=no_permuts,ncol=as.numeric(samplesizes[i,2]))
          for (m in 1:no_permuts) {
            temp_matrix[m,] <- sample(temp_array,as.numeric(samplesizes[i,2]),replace=FALSE)
          }
          temp_matrix <- rbind(rep(input[1,j],as.numeric(samplesizes[i,2])),temp_matrix)
       } else {
         haps <- input[((which(as.numeric(input[2:(dim(input)[1]),j])>=as.numeric(samplesizes[i,2])))+1),1]
         temp_matrix <- matrix(rep(haps,no_permuts),nrow=no_permuts,ncol=length(haps),byrow=TRUE)
         temp_matrix <- rbind(rep(input[1,j],as.numeric(samplesizes[i,2])),temp_matrix)
       }
       record_matrix <- cbind(record_matrix,temp_matrix)  
    }
    write.table(record_matrix,paste(samplesizes[i,1],"_",samplesizes[i,2],"_full_hap_record.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  
  
}


}
