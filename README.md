# duplicating_sequences
Some analyses, like skyline plots through BEAST, need each sequence represented (even duplicates). This code allows you to go from unique sequences to each sequence being represented.

#To use:
Copy the duplicating_sequencing.R function into R and then do the following:

Call the program by duplicating_sequences(working_dir,fasta_file_name,freq_file_name), where:

working_dir == pathway to the folder with your fasta and frequency files e.g. "C:/blahblahblah" 

fasta_file_name == the name of your fasta file e.g. "data.fasta"

freq_file_name == the name of your tab-delimited frequency file, with the sample names (including the > fasta prefix) in the first column, and frequency in the second.

#Example of input:
duplicating_sequences("C:/Users/Folder/","All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta","frequency.txt")

#Example of frequency file:
```
>mtGen01	29
>mtGen02	8
>mtGen03	2
>mtGen04	15
>mtGen05	4
>mtGen13	5
```

#Citation info
This pipeline wouldn't be possible without:

R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/
