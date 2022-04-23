The following step-by-step tutorial shows you how to run a pipeline to liftOver CSREP output (the summary chromatin state map and the probabilistic estimates of state assignments) from hg19 to hg38. We note that this pipeline can only be run and understood after you have followed the main tutorial for CSREP. In the main tutorial, we already showed you how to get the CSREP summary chromatin state maps forwo groups ESC and Brain, in hg19 (note, as a toy example, the data was restricted to chr22). Now, we will show you to to use a separate pipeline to liftOver the output from previous tutorial from hg19 to hg38. This can be useful when you don't want to run the whole CSREP pipeline again for liftOver input data (hg38). 

Here's how the pipeline work, with variable names as how they appear in the ```Snakefile``` and ```config/config.yaml```
- Given a file showing length of chromosomes, we write a bed file with 4 columns: chrom, start, end, chrom_start_end. This file will be given the variable name ```org_segment_fn```. How the file looks:
```
chr1    0       200     chr1_0_200
chr1    200     400     chr1_200_400
chr1    400     600     chr1_400_600
```
- Use ```liftOver``` inside ```liftOver_dir```
