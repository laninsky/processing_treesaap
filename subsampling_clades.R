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
}
