# subsampling_clades
Biogeographic inference can be heavily biased by sampling. This code downsamples by geographic region to try and account for this.

As an input this code takes a file with the clades identifiers in the first column, followed by counts of samples found within those clades for each of the geographic regions in your dataset e.g.
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
PA1	18	1	0	58	1
```

After pasting the entire R function into R (or sourcing it), call it by:
```
subsampling_clades(working_dir,file,no_permuts)
```
Where working_dir is the location with your input file/where output files will be written out, your file name (set up as described above), and the number of permutations you'd like to perform (resampling from the regions with large sample sizes).
