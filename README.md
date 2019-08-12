# haplotyping_sequences
This R-script eats a fasta file, and spits out a list of sequences with identical haplotypes, as well as a fasta file with the unique sequences. It ignores differences between sequences involving IUPAC ambiguity codes, but does take into account indels (i.e. "-"). The R-script can handle sequence breaking over multiple lines.

### Example usage
haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil", "output.fasta")

haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil/test", "All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta")

### Version history
This script was written for:  
Morin, P.A., Foote, A.D., Baker, C.S., Hancock‐Hanser, B.L., Kaschner, K., Mate, B.R., Mesnick, S.L., Pease, V.L., Rosel, P.E. and Alexander, A., 2018. Demography or selection on linked cultural traits or genes? Investigating the driver of low mtDNA diversity in the sperm whale using complementary mitochondrial and nuclear genome analyses. Molecular ecology, 27(11), pp.2604-2619.

This pipeline wouldn't be possible without:

R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

Wickham, H., stringr: Simple, Consistent Wrappers for Common String Operations
