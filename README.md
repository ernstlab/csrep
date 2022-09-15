Hello! Wecome to CSREP- a method that probabilisitcally estimate the chromatin state annotation for a group of related samples. A direct application of CSREP is to calculate the differential chromatin state mappings between two groups with multiple samples. Our manuscript will be made public shortly. 
# Data Availability
## CSREP scores for cell groups from ROADMAP
The data of summary chromatin state maps outputted by CSREP for all cell groups in ROADMAP and Epimap is available at https://public.hoffman2.idre.ucla.edu/ernst/2K9RS/csrep/
Below is the file structure of this download folder
```
|__ roadmap
|__|__ <cell_group>
|__|__|__ hg19
|__|__|__|__ summary_state_track.bed.gz: This file shows the summary chromatin state map for the group. The file is designed such that it had a header that can be read into UCSC genome browser
|__|__|__|__ state_assign_matrix: inside this folder, you will see the data of probabilistic estimates of state assignments for each genomic position in the genome. There are 316 files, each  corresponding to a window of 10Mb (or less, if the file belongs to the end of a chromosome). A file chr1_0_avg_pred.txt.gz represents region chr1: 0: 9,999,999 bp. File chr5_15_avg_pred.txt.gz represents region chr5: 150Mp- 150Mb+9,999,999bp. Each line in the file represents a 200-bp window in the genome. If it's too confusing, just have a look at the file and you will figure out. 
|__|__|__|__ ucsc_tracks: inside this folder, we stored the bigwig and bigbed files used for visualizing the summary state maps resulted from CSREP on UCSC Genome Browser
|__|__|__ hg38
|__|__|__|__ summary_state_track.bed.gz: This file is lifted-over from the corresponding output in hg19
|__|__|__|__ state_assign_matrix: inside this folder, you will see 23 files corresponding to 23 chromosomes (1-22,X, we do not calculate output for chrY). Data here is actually lifted over from the CSREP's output in hg19. Implementations and details about how we get CSREP output by lifting over the output from another assembly can be found at <hahahaha fill in>  
|__|__|__|__ ucsc_tracks: inside this folder, we stored the bigwig and bigbed files used for visualizing the summary state maps resulted from CSREP on UCSC Genome Browser
|__ epimap
|__|__ <cell_group>
|__|__|__ hg19
|__|__|__|__ summary_state_track.bed.gz: Same structure as in Roadmap (mentioned above)
|__|__|__|__ state_assign_matrix: Same structure as in Roadmap (mentioned above)
|__|__|__|__ ucsc_tracks: Same structure as in Roadmap (mentioned above)
|__|__|__ hg38: 
|__|__|__|__ summary_state_track.bed.gz: Same structure as in Roadmap (mentioned above)
|__|__|__|__ state_assign_matrix: Same structure as in Roadmap (mentioned above)
|__|__|__|__ ucsc_tracks: Same structure as in Roadmap (mentioned above)

```

## Viewing summary chromatin state maps on USCS Genome Browser
Users can easily view the provided summary chromatin state maps (state assignment probabilities and state annotations) for cell groups in Roadmap and Epimap on USCS Genome Browser, using the link to the track hub: https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//csrep/hub.txt


## Other data
For any other data related to the paper, please email Prof. Jason Ernst or grad student Ha Vu. 

# Note from the authors:
If you run into any problems running CSREP or following the tutorial accompanying the software, please contact grad student Ha Vu via email havu73@ucla.edu. We will be happy to assist and hear your feedback to make the software more fool-proof and user-friendly. 

# Installing CSREP
## Software requirement
All our package specifications needed for running CSREP are provided in file ```env/csrep_env.yml``` (Note: path to this file is relative to the current directory that contains this README file).

 
## How to install CSREP
- Clone this github repository
- Replicate our conda environment: ```conda env create -f env/csrep_env.yml``` (Note: path to this file is relative to the current directory that contains this README file)

## Test data
To show you how CSREP can be run, we provided an example that include 18-chromatin-state segmentation data for 5 samples of ESC groups and 7 samples of Brain group from <a href="https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state">Roadmap Epigenomics Project</a>. We restricted the example input to include only data from chromosome 22.

# Run CSREP 
## Snakemake 
We highly recommend using Snakemake for running CSREP. The conda environment provided in ```env/csrep_env.yml``` already includes snakemake as a requriement. We have written the Snakefile code for which you can easily modify parameters in file ```config/config.yaml``` in order to apply CSREP for your own purposes. Reasons that we highly recommend using snakemake include: (1) If you have computing cluster available, you can easily train multiple models for CSREP in parallel through different jobs using Snakemake --> huge reduction in runtime. Snakemake gracefully manages the parallel jobs and automatically rerun jobs that were left unfinished. Additionally, Snakemake manages pipeline in such a way that it can run jobs from where the pipeline were unfinished or changed, therefore saving users lots of time babysitting jobs. (2) Snakemake ensures the consistency in computing environment, therefore ensuring the reproducibility of our results. 

In order to run CSREP either to **summarize** the chromatin state maps for one or multiple group(s) of samples, or to **calculate the differential** chromatin scores between two groups, you can follow two steps: 
(1) prepare input data files/folders
(2) modify the config file that we provided with user-input parameters (this step is very quick once you understand the meaning of different parameters in ```config/config.yaml```). 
Details about how to run CSREP are below:

## Step 1: Prepare input data
Users need to prepare the following files/folders:
- **Chromatin state segmentation annotation files** for samples associated with:
(1) the groups that we want to calculate the summary chromatin state maps (CSREP can run for multiple potentially overlapping groups concurrently) 
or with 
(2) each of the two groups that we want to calculate differential scores. 
Chromatin state segmentation files should be in BED format, with each chromatin state defined at 200bp resolution (default settings of ChromHMM or SegWay). Output files from ChromHMM are supported as input to CSREP. We provided example input files in ```testdata``` folder, and they will also be used in the tutorial.
- **A text file listing out the samples associated with each group we would like to run CSREP on**: For each group of samples, you will need to prepare a file where each line shows a sampleID. For example, the ESC group from Roadmap has 5 samples: ```E003, E008, E014, E015, E016```, each sample has a segmentation file ```<sampleID>_chr22_core_K27ac_segments.bed.gz``` as input data in our tutorial. Then, the sample file for ESC look like: 
```
E003
E008
E014
E015
E016
```
In step 2 below, we will note where this file should be placed. For now, users can just focus on creating this file. Sorry for the inconvenience. 
- **A BED file specifying the length of each chromosome in the genome of interest** (hg19, hg38, mm10, etc.). The first column shows the chromosome, the second column is 0 in all rows, the third column shows the length of the chromosome. We provided code ```utils/get_chrom_length_from_segmentation.sh``` to create this file of chromosome length based on the input segmentation data of one sample. Details about how to run this code will be provided in the ```utils/README.md``` file. 

## Step 2: Modify parameters to CSREP
CSREP is designed such that you can simply modify the ```config/config.yaml``` file with desired parameter values, then run snakemake to get results from CSREP. CSREP allows users to calculate the summary chromatin state maps for multiple groups of samples simultaneously.

Below, we present the parameters and the file structures outputted by CSREP. We also provide a tutorial for running CSREP on our test data. If reading the following list of parameters gets too confusing, you can skim through it, then go to the tutorial and look back at this list for reference. We also provide file ```cheatsheet_csrep_manual.pdf``` to help you better understand the parameters. 

Parameters in ```config/config.yaml``` include:
- ```is_calculate_diff_two_groups```: **0** if you want to calculate the summary (representative) chromatin state map for >=1 groups of samples. **1** if you want to calculate the differential chromatin scores between two groups of samples.
- ```all_cg_out_dir```:  The folder where the **output** data (representative/differenntial chromatin state maps) for all groups of samples are stored. Within this folder, each subfolder will correspond to a group (for the summary task) or to a pair of groups (for the differential state calculation task). You decide through ```all_cg_out_dir``` where this data is stored.
- ```raw_user_input_dir```: the folder where all **input segmentation data files** of all samples, regardless of their group memberships, are stored.
- ```all_ct_segment_folder```: A folder where we store the segmentation data of all samples, after processing from the ```raw_user_input_dir```. The tutorial will be helpful in understanding the role of this folder, if you care. You decide through ```all_ct_segment_folder``` where this data is stored.
- ```raw_segment_suffix```: The suffix of input samples' segmentation files. In the tutorial, inside input folder ```test_data/raw_data``` (corresponding to ```raw_user_input_dir```), each sample's segmentation file name is of format ```<sampleID>_chr22_core_K27ac_segments.bed.gz```, then ```raw_segment_suffix``` would be ```_chr22_core_K27ac_segments.bed.gz```.
- ```training_data_folder```: the folder where CSREP stores chromatin state segmentation data used for training, which corresponds to 10% of the genome, for each sample. You decide through ```training_data_folder``` where this data is stored.
- ```chrom_length_fn```: path to the BED file that stores the length of chromosomes. This is useful in sampling the training regions for CSREP. You can use ```utils/get_chrom_length_from_segmentation.sh``` to generate this file, or you can write it yourself. 
- ```sample_genome_fn```: a BED file where the locations of training regions for CSREP are stored. You decide through ```sample_genome_fn``` where this data is stored.
- ```chromhmm_state_num```: number of chromatin states that are in the model for each sample. 
- ```train_mode_list```: the list of training mode that you would like results for. ```multi_logistic``` for CSREP and ```baseline``` for base_count, as presented in the paper. 
- ```cell_group_list```: the list of group names of groups of samples in that we would like to calculate the representative/differential chromatin state maps for. If ```is_calculate_diff_two_groups``` is set to 1 (meaning you want to calculate differential chromatin scores), then ```cell_group_list``` should specify two group names, otherwise the program may inform you of an error.

# Tutorial
Please see this <a href="https://github.com/ernstlab/csrep/blob/master/tutorial.md">link</a>. If the info in section 'How to run CSREP' confuse you, you can try following our tutorial, it will make things a lot easier to understand.

# liftOver results to another assembly
Given the output of summary/differential chromatin state maps from CSREP in one assembly, you can liftOver such results into another assembly. We provide source code and detailed readme of a separate pipeline for this. The pipeline is designed such that we will first lift the segmentation one bin at at time, and then eradicate regions in the destination assembly that get mapped from multiple regions in the original assembly. We used this pipeline to create the summary chromatin state maps for 11 cell groups from Roadmap from hg19 (produced by CSRE) to hg38. The pipeline can be found <a href="https://github.com/ernstlab/csrep/tree/master/liftOver_csrep_output">here</a>. The accompanying readme/tutorial is <a href="https://github.com/ernstlab/csrep/blob/master/liftOver_csrep_output/README.md">here</a>

# License
All code is provided under the MIT Open Acess License Copyright 2021 Ha Vu and Jason Ernst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.