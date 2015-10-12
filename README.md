# haplotyping_sequences
This R-script eats a fasta file, and spits out a list of sequences with identical haplotypes, as well as a fasta file with the unique sequences. It ignores differences between sequences involving IUPAC ambiguity codes, but does take into account indels (i.e. "-"). The R-script can handle sequence breaking over multiple lines.

#Example usage
haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil", "output.fasta")

haplotyping_sequences("C:/Users/Alana/Dropbox/Mitogenome_Phil/test", "All_Pmac_Unique_Haps_13gene_Concat_HapAssign.fasta")
