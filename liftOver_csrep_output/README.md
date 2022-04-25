The following step-by-step tutorial shows you how to run a pipeline to liftOver CSREP output (the summary chromatin state map and the probabilistic estimates of state assignments) from hg19 to hg38. We note that this pipeline can only be run and understood after you have followed the main tutorial for CSREP. In the main tutorial, we already showed you how to get the CSREP summary chromatin state maps forwo groups ESC and Brain, in hg19 (note, as a toy example, the data was restricted to chr22). Now, we will show you to to use a separate pipeline to liftOver the output from previous tutorial from hg19 to hg38. This can be useful when you don't want to run the whole CSREP pipeline again for liftOver input data (hg38). 

Here's how the pipeline work, with variable names as how they appear in the ```Snakefile``` and ```config/config.yaml```
- Assume we want to convert the output from CSREP from ```org_assembly``` to ```end_assembly```. Since the the liftOver coordinate file (downloaded from UCSC genome browser, stored at ```liftOver```) is named ```hg19ToHg38.over.chain.gz```, in the configuration (```config/config.yaml```), we specified ```org_assembly``` as ```hg19``` and ```end_assembly``` as ```Hg38```. Note the capitalizations of ```hg19``` and ```Hg38``` to match with the file name ```hg19ToHg38.over.chain.gz```. 
- Given a file showing length of chromosomes, we write a bed file with 4 columns: chrom, start, end, chrom_start_end. This file will be given the variable name ```org_segment_fn```. How the file looks:
```
chr1    0       200     chr1_0_200
chr1    200     400     chr1_200_400
chr1    400     600     chr1_400_600
```
- We will specify ```end_segment_dir``` in the configuration as a folder where the pipeline will store the files that show 1-1 mappings regions between the two assemblies. If you are confused about what to put it, just keep it as is how we put it in our current file ```config/config.yaml```. 
- Use ```liftOver``` inside ```liftOver_dir``` to convert ```org_segment_fn``` from ```org_assembly``` to ```end_assembly```. Then, we will sort the lifted-over files by genomic coordinates, and get rid of any regions in the ```end_assembly``` that were mapped from multiple different regions from ```org_assembly```. This process will produce many temporary files inside this current working directory, please let snakemake does its jobs and it will clean up all the temp. files later. In the end, file ```destOrg_fn``` will be produced, with columns: chrom, start, end, orgC_orgS_orgE. The first 3 columns show the genomic coordinates in ```end_assembly```. The last column shows the 1-1 mapped chrom_start_end coordinate of the region in ```org_assembly```. How the file ```destOrg_fn``` looks:
```
chr1    10000   10200   chr1_10000_10200
chr1    10200   10400   chr1_10200_10400
chr1    10400   10600   chr1_10400_10600
```
This means, for example, regions chr1:10,000-10,200 in ```org_assembly```  (last column) is mapped to chr1:10,000-10,200 in ```end_assembly``` (first 3 columns)

- From ```destOrg_fn```, we will then produce the file ```orgDest_fn```, which is the inverse mapping bed file of ```destOrg_fn```. The file will have 4 columns: chrom, start, end, destC_destS_destE. The first 3 columns show the genomic coordinates in ```org_assembly```. The last column shows the 1-1 mapped chrom_start_end coordinate of the region in ```end_assembly```. How the file ```orgDest_fn``` looks:
```
chr1    10000   10200   chr1_10000_10200
chr1    10200   10400   chr1_10200_10400
chr1    10400   10600   chr1_10400_10600
```
This means, for example, regions chr1:10,000-10,200 in ```org_assembly```  (first 3 columns) is mapped to chr1:10,000-10,200 in ```end_assembly``` (last column). 

- Now, we will try to convert the CSREP output. The procedure simply takes in the CSREP summary data, and the 1-1 mapping of genomic regions that we produced in ```orgDest_fn```, and will map genomic bins from CSREP output (in ```org_assembly```) to genomic bins in ```dest_assembly```. 
To run this conversion code, we will need to specify some varibles in the configurations (```config/config.yaml```). Some variables in this pipeline will correspond to those in the CSREP pipeline (refer to the main readme and the CSREP tutorial as you go along this tutorial, if you are confused):
(1) We specified ```all_cg_out_dir``` in CSREP pipeline, where there are subfolders corresponding to different groups of samples. The exact file path should be copied to the config file of this liftOver pipeline to the variable name ```org_all_ct_folder```. This variable refers to the folder from ```org_assembly``` where we can find subfolders of the cell groups that we obtained CSREP output for. 
(2) We specified ```cell_group_list``` in CSREP pipeline, specifying the list of cell groups that we would like to obtain CSREP summary/differential data for. In this liftOver pipeline, we will also fill in the variable ```cell_group_list``` in the configuration. This will specify the list of cell groups that we would like to liftOver CSREP data for, therefore, ```cell_group_list``` should be a subset of ```cell_group_list```. 
(3) We specified ```train_mode_list``` in CSREP pipeline, which can be a subset of ```['multi_logistic', 'baseline']```. In this liftOver pipeline, we will also fill in the varible ```train_mode_list``` based on what methods' summary output (CSREP or baseline or both) that you would like to liftOver. Yes, we can liftOver the output of base_count (as named in the manuscript) too. 
(4) We specified ```chromhmm_state_num``` as the number of chromatin states in the input annotations. We will do the same here. 
(5) In this liftOver pipeline, we will specify ```gene_reg_list``` as the list of genomic regions in the genome. We set them to a list of two regions now in the configuration, you can keep it as is.
(6) The output will include: (1) the summary chromatin state map in ```dest_assembly```, and the state assignment matrix in ```dest_assembly```, all will be stored inside folder ```dest_all_ct_folder``` (users, specify this folder in the ```config/config.yaml```). The structure of the output are as follows:
```dest_all_ct_folder```
|__ ```<cell_group>```
|__|__ ```<train_mode>```: multi_logistic (CSREP) or baseline (base_count)
|__|__|__ representative_data
|__|__|__|__ summary_state_track.bed.gz: file of the summary chromatin state map for the ```<cell_group>```
|__|__|__|__ state_assign_matrix: folder where the data of probabilities of state assignments across the genome in ```dest_assembly``` are stored
|__|__|__|__|__ ```<chrom>_liftOver_probState.txt.gz```: Each file shows the data of state assignment matrix for one chromosome. 

