Hello! Wecome to CSREP- a method that probabilisitcally estimate the chromatin state annotation for a group of related samples. A direct application of CSREP is to calculate the differential chromatin state mappings between two groups with multiple patterns. Our manuscript will be made public shortly. 
# Installing CSREP
## Software requirement
In order for CSREP to work, we need:
- bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
We create an environment that is compatible with our program, which you install as a conda environment: ```conda env create -f env/csrep_env.yml```

## Test data
To show you how CSREP can be run, we provided an example that include 18-chromatin-state segmentation data for 5 samples of ESC groups and 7 samples of Brain group from Roadmap Epigenomics Project (https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state). We restricted the example input to include only data from chromosome 22.
 
## How to install CSREP
- Clone this github repository
- Replicate our conda environment: ```conda env create -f env/csrep_env.yml``` (Note: path to this file is relative to the current directory that contains this README file)
- Install bedtools. 

# Run CSREP 
## Snakemake 
We highly recommend using Snakemake for running CSREP. The conda environment provided in ```env/csrep_env.yml``` is already compatible to apply snakemake for running CSREP. We have written the Snakefile code for which you can easily use for your own purposes by modifying parameter in file ```config/config.yaml```. Reasons that we highly recommend using snakemake include: (1) If you have computing cluster available, you can easily train multiple models for CSREP in parallel through different jobs using Snakemake --> huge reduction in runtime. Snakemake gracefully manages the parallel jobs and automatically rerun jobs that were left unfinished. Additionally, Snakemake manages pipeline in such a way that it can run jobs from where the pipeline were unfinished or changed, therefore saving users lots of time babysitting jobs. (2) Snakemake ensures the consistency in computing environment, therefore ensuring the reproducibility of our results. 

In order to run CSREP either to summarize the chromatin state maps for a group of sample or to calculate the differential chromatin landscape between two groups with multiple samples, you can easily modify the config file that we provided with parameters. Details about how to run sciddo are below

## Input data
In order to run CSREP, users need to prepare the following files:
- Chromatin state segmentation annotation files for each of the samples involved in (1) the groups that we want to calculate the summary chromatin state maps (CSREP can run for multiple groups concurrently) or in (2) each of two groups (for the calculation of differential scores task). CSREP supports text-based input and output data formats. Chromatin state segmentation files in BED format are support at 200bp resolution, which is the resolution used in most of chromatin state annotation produced by ChromHMM or Segway. Output files from ChromHMM are supported as input to CSREP. We provided example input files which will be used in the attached tutorial.
- A text file listing out the samples associated with each group of samples we would like to run CSREP on. For each group of samples, you will need to prepare a text file with each line corresponding to a sample in such group. For example, the ESC group from Roadmap has 5 samples: ```E003, E008, E014, E015, E016```, each sample has a segmentation file ```<sample_code>_chr22_core_K27ac_segments.bed.gz``` as input data in our tutorial. Then, the sample file for ESC look like: 
```
E003
E008
E014
E015
E016
```
- A file in BED format specifying the length of each chromosome in the genome of interest (hg19, hg38, mm10, etc.). The first column shows the chromosome, the second column is 0 in all rows, the third column shows the length of the chromosome. We provided code ```utils/get_chrom_length_from_segmentation.sh``` to creat this file of chromosome length based on the input segmentation data of one sample. Details about how to run this code will be provided in the ```utils/README.md``` file. 

## Modify parameters to CSREP
Right now, we create CSREP such that you can simply modify the ```config/config.yaml``` file with your specified parameters, then run snakemake to get results from CSREP. CSREP is designed such that you can calculate the summary chromatin state map for multiple groups of samples at the same time. Below, we present the parameters and the file structures outputted by CSREP, please also refer to the file ```csrep_file_structure.pdf``` to help you better understand the parameters. We also provided a tutorial to help you get started on CSREP and more easily run CSREP for your purposes. Parameters in ```config/config.yaml``` include:
- ```is_calculate_diff_two_groups```: 0 if you want to calculate the summary (representative) chromatin state map for >=1 groups of samples. 1 if you want to calculate the differential chromatin scores between two groups of samples.
- ```all_cg_out_dir```:  The output folder where the output data (representative/differenntial chromatin state maps) for all groups of samples are stored. Within this folder, each subfolder will correspond to a group (for the summary task) or to a pair of groups (for the differential state calculation task). You decide through ```all_cg_out_dir``` where this data is stored.
- ```raw_user_input_dir```: the folder where all input segmentatil data files of all samples, regardless of their group memberships, are stored.
- ```all_ct_segment_folder```: A folder where we store the segmentation data of all samples, after processing from the ```raw_user_input_dir```. The tutorial will be helpful in understanding the role of this folder, if you care. You decide through ```all_ct_segment_folder``` where this data is stored.
- ```raw_segment_suffix```: The suffix of input samples' segmentation files. In the tutorial, inside input folder ```test_data/raw_data``` (corresponding to ```raw_user_input_dir```), each samples' segmentation file name is of format ```<sample_ids>_chr22_core_K27ac_segments.bed.gz```, then ```raw_segment_suffix``` would be ```_chr22_core_K27ac_segments.bed.gz```.
- ```training_data_folder```: the folder where we will store chromatin state segmentation data used for training, which corresponds to 10% of the genome, for each sample. You decide through ```training_data_folder``` where this data is stored.
- ```chrom_length_fn```: path to the BED file that stores the length of chromosomes. This is useful in sampling the training regions for CSREP. You can use ```utils/get_chrom_length_from_segmentation.sh``` to generate this file, or you can write it yourself. 
- ```sample_genome_fn```: a BED file where the locations of training regions for CSREP are stored. You decide through ```sample_genome_fn``` where this data is stored.s
- ```chromhmm_state_num```: number of chromatin states that are in the model for each sample. 
- ```train_mode_list```: the list of training mode that you would like results for. ```multi_logistic``` for CSREP and ```baseline``` for base_count, as presented in the paper. 
- ```cell_group_list```: the list of group names of groups of samples in that we would like to calculate the representative/differential chromatin state maps for. If ```is_calculate_diff_two_groups``` is set to 1 (meaning you want to calculate differential chromatin scores), then ```cell_group_list``` should specify two group names, otherwise the program may crash

# Data Availability
## Male-Female comparision
The data of histogram for differential scores for Male-minus-Female is available at https://public.hoffman2.idre.ucla.edu/ernst/2K9RS/Male_minus_Female/multi_logistic/histogram/ 
## CSREP scores for cell groups from ROADMAP
The data of genome-wide CSREP score for all cell groups in ROADMAP is available at https://public.hoffman2.idre.ucla.edu/ernst/2K9RS/
The 11 cell groups's subfolders include: blood Blood_other Digestive Brain ESC ES-deriv Heart iPSC Muscle Skin Sm_Muscle. To access the csrep score for a group, enter `<cell_group_name>`/`csrep`
## Other data
For any other data related to the paper, please email Prof. Jason Ernst or grad student Ha Vu. 
