# subsampling_clades v.1.0.0
Originally at: https://github.com/laninsky/subsampling_clades

Biogeographic inference can be heavily biased by sampling. This code downsamples by geographic region to try and account for this.

As an input this code takes a file with the clades identifiers in the first column (whitespace delimited), followed by counts of samples found within those clades for each of the geographic regions in your dataset e.g.
```
Clades	Atlantic	GoM	Mediterranean	Pacific	Indian
A1	11	4	0	0	0
A2	6	6	0	0	0
M1	0	0	4	0	0
P1	0	0	0	2	0
P2	0	0	0	3	0
P3	0	0	0	8	0
P4	0	0	0	8	0
P5	0	0	0	28	0
P6	0	0	0	10	0
P7	0	0	0	7	0
PA1	8	0	0	0	0
PAI	10	1	0	58	1
```

After pasting the entire R function into R (or sourcing it), call it by:
```
subsampling_clades(working_dir,file,no_permuts)
```
Where working_dir is the location with your input file/where output files will be written out, your file name (set up as described above), and the number of permutations you'd like to perform (resampling from the regions with large sample sizes).

The program ticks through each region's sample size (from small to big), and subsamples any region which has a sample size larger than the current region. For the example above, the most extreme case is the Indian Ocean where only 1 individual was sampled. Using a sample size of one, the program resampled each other ocean. The raw output (labelled as Region_samplesize_full_clade_record.txt e.g. Indian_1_full_clade_record.txt has iterations in rows and the clades taken from oceans in columns with ocean headers in the top row e.g. (first four iterations)
```
Atlantic GoM Mediterranean Pacific Indian
A2 PAI M1 P3 PAI
PAI A2 M1 PAI PAI
PA1 A1 M1 P7 PAI
PAI PAI M1 P5 PAI
```
The program also summarizes the proportion of iterations where each clade was present for each ocean (labelled as Region_samplesize_summary_record.txt e.g. Indian_1_summary_record.txt ) e.g.
```
Clades Atlantic GoM Mediterranean Pacific Indian
A1 0.293 0.378 0 0 0
A2 0.167 0.534 0 0 0
M1 0 0 1 0 0
P1 0 0 0 0.012 0
P2 0 0 0 0.017 0
P3 0 0 0 0.079 0
P4 0 0 0 0.068 0
P5 0 0 0 0.225 0
P6 0 0 0 0.085 0
P7 0 0 0 0.057 0
PA1 0.235 0 0 0 0
PAI 0.305 0.088 0 0.457 1
```
e.g. only 29.3% of the time when pulling just a single sample from the Atlantic, does this correspond to clade A1.

If the sample size of a region is smaller or equal to the current sample size being used to subsample other regions, then the observed clades are output e.g. for the Gulf of Mexico sample size (11), GoM_11_full_clade_record.txt:
```
Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic Atlantic GoM GoM GoM Mediterranean Pacific Pacific Pacific Pacific Pacific Pacific Pacific Pacific Pacific Pacific Pacific Indian
PAI PA1 PAI PA1 PA1 PAI A1 A2 A1 A1 A1 A1 A2 PAI M1 PAI P3 PAI P4 P5 PAI P3 PAI PAI P5 PAI PAI
PA1 A2 PAI A2 PAI PA1 A1 PAI PA1 A1 PAI A1 A2 PAI M1 PAI P5 PAI P5 PAI PAI P5 PAI PAI PAI PAI PAI
PA1 PAI A1 A1 PAI A1 PAI PA1 A1 A1 A1 A1 A2 PAI M1 P2 P4 P7 P5 P3 PAI P5 PAI P6 P2 P4 PAI
A1 PA1 PA1 A1 A1 A1 A1 PAI PAI A2 PA1 A1 A2 PAI M1 P5 P1 P4 PAI PAI P6 P5 P5 PAI PAI PAI PAI
```
Clades are not downsampled for GoM, Med or Indian. The clades output are the same found in the observed input file (A1, A2 and PAI for GoM, M1 for Med, PAI for Indian).

Example of summary file (GoM_11_summary_record.txt) for GoM's sample size of 11 (used to resample Atlantic and Pacific which have larger sample sizes than GoM):
```
Clades Atlantic GoM Mediterranean Pacific Indian
A1 0.989 1 0 0 0
A2 0.912 1 0 0 0
M1 0 0 1 0 0
P1 0 0 0 0.167 0
P2 0 0 0 0.243 0
P3 0 0 0 0.523 0
P4 0 0 0 0.533 0
P5 0 0 0 0.947 0
P6 0 0 0 0.625 0
P7 0 0 0 0.477 0
PA1 0.971 0 0 0 0
PAI 0.99 1 0 1 1
```
The PAI clade is found in both the Atlantic and Pacific in the vast majority of the permutations, despite only 11 samples being re-sampled from each ocean in each permutation. The fact that we observe this clade for both of these oceans is probably not solely due to the fact that we have larger sample sizes from these oceans. Some of the other clades (e.g. P1 in the Pacific) could potentially have been missed in other oceans due to them having smaller sample sizes, because we only observed this clade in roughly 17% of simulations when we downsample to 11 individuals in the Pacific.

# Version history
v1.0.0 ready to go. This version first published with:  
Morin, P.A., Foote, A.D., Baker, C.S., Hancock‐Hanser, B.L., Kaschner, K., Mate, B.R., Mesnick, S.L., Pease, V.L., Rosel, P.E. and Alexander, A., 2018. Demography or selection on linked cultural traits or genes? Investigating the driver of low mtDNA diversity in the sperm whale using complementary mitochondrial and nuclear genome analyses. Molecular ecology, 27(11), pp.2604-2619.

If you use this script please cite the publication above and the following:  
Alexander, A. 2018. subsampling_clades v.1.0.0. Available from: https://github.com/laninsky/subsampling_clades

This script also wouldn't be possible without:  
R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/
