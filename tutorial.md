The following step-by-step tutorial shows you how to run CSREP to either (1) calculate the summary chromatin state map for a group of samples, or (2) calculate the differential chromatin state map between two groups of samples. Right now, we implement CSREP such that users only need to modify paramters in a config file ```config/config.yaml```
We provided some test data for ease of following this tutorial (in folder ./testdata/raw_data). These are data of chromatin state segmentation in chromosome 22 for 5 samples of the ESC group and 7 samples of the Brain group from Roadmap Epigenetic Project. We chose to include only chromosome 22 to make sure our example can be run in reasonable time for demonstration purposes. We used the segmentation data from the [18-state model](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state), which was trained on 6 histone modification [marks: H3K4me3, H3K27ac, H3K4me3, H3K36me3, H3K27me3, H3K9me3](https://egg2.wustl.edu/roadmap/figures/extendedData/Figure_ED2.jpg) 
# File structure for reference in running CSREP 
**You can briefly read through sections of step 0 and 1 first before going back to reading this section in details. Also, this can serve as an appendix in your first exploration run of CSREP with snakemake**

The following text shows you the file structure in which CSREP organizes input and output files. We write them in the terms of the variable names in the config file ```config/config.yml```. You can look at this file structure for reference, which will be helpful when you fill in parameters for CSREP in the config file

```testdata/raw_data``` corresponding to ```raw_user_input_dir```: the folder where all input segmentation data files of all samples, regardless of their group memberships, are stored. Inside this folder: 
\|\_\_```<sample_id>```: each of these folders correspond to a sample. In the example, there are 12 of these folders. Inside each folder: 
\|\_\_\|\_\_```<sampleID><input_filename_suffix>```: This file corresponds to one biosample's input segmentation input. In the example provided with the tutorial, each of these files are in the format ```<sampleID>_chr22_core_K27ac_segments.bed.gz```, so ```input_filename_suffix``` will then be ```_chr22_core_K27ac_segments.bed.gz```. **Users will prepare data in this folder. The format of this file requires, at the minimum, the first 4 columns to be of the form:**
```
<chrom>	<start_bp>	<end_bp>	<state>
```
And example of the file format is as follows:
```
chr10   0       94400   E18
chr10   94400   112000  E13
chr10   112000  117200  E18
``` 
**Note about state format**: The states are usually of the form ```E<state_index>``` (where ```<state_index>``` is one-based), which is the default format if the segmentation is being done by ChromHMM. You can also have your states being of the form ```<one_character><state_index>```, examples for state ```1```: ```E1, U1, S1```, etc. Or, sometimes your states may be the full state names or the state mnenomic (ex: ```TSSA, Quies, EnhA```, etc.) --> if this is the case, your ```state_annot_fn``` should have at least two columns: ```state``` and ```mnenomic``` to specify the meaning of states and their indices. The state format is important because CSREP will try to convert the state names/annotations in the input into a state_number system. If this step is wrong, it may corrupt the entire CSREP pipeline. Details of the ```state_annot_fn``` file are provided below. 

```testdata/state_annot.txt``` corresponding to ```state_annot_fn```. This file specifies different chracteristics of the states. **Users will prepare data in this File. The format of this file requires, at the minimum, columns ```state```, ```itemRgb``` and ```mnenomic```**. Column ```state``` should show the states' indices (one-based), column ```mnenomic``` should show the state names (Many repositories such as [Epimap](http://compbio.mit.edu/epimap/) provides raw data of state names in the mnenomic forms such as ```TSSA, Quies, EnhA```, etc.). Column ```itemRgb``` shows the state colors, in ```r,g,b``` format. We provide an example state annotation file ```testdata/state_annot.txt``` for your reference.

```testdata/csrep_output``` corresponding to ```<all_cg_out_dir>```:  The output folder where the output data (representative/differenntial chromatin state maps) for all groups of samples are stored. All data inside this folder will be produced as part of the CSREP pipeline, except for ```sample.list``` files which users will have to prepare. 

```testdata/roadmap_18state_chromlength.bed``` corresponding to ```chrom_length_fn```: bed file showing the length of each chromosome. **This file can be produced in Step 0 by users.**

```testdata/sample_genome.bed.gz``` corresponding to ```sample_genome_fn```: where the data of regions that we sample for training data are stored. Data in this folder is produced by CSREP.

\|\_\_ ```Brain```: the folder where we can get summary chromatin state map data for samples in ```Brain``` group. The structure of this folder is similar to that of ```ESC``` folder as shown belows.\
\|\_\_ ```ESC```: the folder where we can get summary chromatin state map data for samples in ```ESC``` group. \
\|\_\_\|\_\_ ```sample.list```: file storing all the sampleIDs in the ```ESC``` group. Each sampleID is on a separate line (```E003, E008, E014, etc.```). **This is a file where the users have to create and place them within this folder in order for Snakemake to work.**\
\|\_\_\|\_\_ ```CSREP```: where output from CSREP's summary chromatin state maps are stored\
\|\_\_\|\_\_\|\_\_ ```representative_data```\
\|\_\_\|\_\_\|\_\_\|\_\_ ```average_predictions```: folder containing the representative chromatin state maps for the group. This is the main output that you are interested in for summarizing a group of samples' chromatin state maps. Each file in this folder has the format ```chr<chrom>_<region_index, 0-based>_avg_pred.txt.gz```, corresponding to a region of at most 10M bp on the chromosome. Each line in the file corresponds to one bin of the chromatin state map (typically 200bp, therefore, a file containing data for 10Mbp window will have 50K lines, each showing the chroamtin state map for a 200bp region). File ```chr22_0_avg_pred.txt.gz``` shows results in the first 10Mbp window of chromosome 22. File ```chr22_5_avg_pred.txt.gz``` shows results in the last 1304400 bp in the chromosome, since chr22 has length 51304400 bp. The first line of each file shows headers (chromatin states), following lines show the probabilities of state assignment at each genomic bin. One file in this folder will be name ```summary_state_track.bed.gz```, showing the state that have the max assignment probabilities (based on files ```chr<chrom>_<region_index, 0-based>_avg_pred.txt.gz```) at each genomic position. This file can be read into UCSC genome browser if the user-provided  ```state_annot_fn``` are correctly formatted with the required columns (see our provided example at ```testdata/state_annot.txt``` for a sample).  \
\|\_\_\|\_\_\|\_\_\|\_\_ ```pred_<sampleID>```: example: ```pred_E003```, a folder containing the data of chromatin state map prediction of the sample (```E003```) by training a multi logistic regression model using input from other samples in the group (see the paper for more details). These folders are redundant and may not need to be looked at for the final output of representative chromaitn state map, therefore, they can be deleted manually by users. We keep these folders as a design choice because keeping them is beneficial for running snakemake for such a data-intensive task. \
\|\_\_\|\_\_ ```baseline```: where output from the base_count method for summary chromatin state maps are stored, if users specify through ```config/config.yml```  \
\|\_\_\|\_\_\|\_\_ ```representative_data```\
\|\_\_\|\_\_\|\_\_\|\_\_ ```average_predictions```: Formatted similarly as the counterpart in the ```CSREP``` parent folder. The results show calculation of the summary chromatin state maps by calculating the frequencies of each state across samples at each position (Paper/methods/base_count)\
\|\_\_ ```ESC_minus_Brain```: Where the data of differential chromatin state maps between the two groups ```ESC``` and ```Brain``` are stored.\
\|\_\_\|\_\_ ```CSREP```: data of differential chromatin state maps calculated using CSREP. Within this folders, each region of at most 10Mbp are presented in one file. All file names and file format as similar to that of the folder containing representative chromatin state maps for one group fo samples. \
\|\_\_\|\_\_ ```baseline```: data of differential chroamtin state maps calculated using base_count.\
```testdata/training_data``` corresponding to ```training_data_folder```: the folder where we will store chromatin state segmentation data used for training, which corresponds to 10% of the genome, for each sample. Data in this folder is produced by CSREP.\




# Step 0: Preparing some input
## Getting data of chromosome length
One input to the CSREP framework is a bed file that show the length of each chromosome. An example file is is as follows: 
```
chr1    0       249250600
chr2    0       243199200
chr3    0       198022400
chr4    0       191154200
chr5    0       180915200
chr6    0       171115000
chr7    0       159138600
chr8    0       146364000
chr9    0       141213400
chr10   0       135534600
chr11   0       135006400
chr12   0       133851800
chr13   0       115169800
chr14   0       107349400
chr15   0       102531200
chr16   0       903546000
chr17   0       811952000
chr18   0       780772000
chr19   0       591288000
chr20   0       630254000
chr21   0       481298000
chr22   0       513044000
chrX    0       155270400
chrY    0       593734000
```
Note: We want the length of the chromosome is a multiple of 200, which is the resolution of the chromatin state segmentation. We recommend using our code ```utils/get_chrom_length_from_segmentation.sh``` to produce this file. The command-line format for running this file is: 
```
./utils/get_chrom_length_from_segmentation.sh <input segmentation file> <output file>
```
For our example, you can run: 
```
./utils/get_chrom_length_from_segmentation.sh ./testdata/raw_data/E072_chr22_core_K27ac_segments.bed.gz ./testdata/roadmap_18state_chromlength.bed
```

## Prepare the list of samples for each group
Inside ```./testdata/raw_data```, there are 12 segmentation files in the form ```<sampleID><input_filename_suffix>```, for example: 
```
E003_chr22_core_K27ac_segments.bed.gz
E008_chr22_core_K27ac_segments.bed.gz
E014_chr22_core_K27ac_segments.bed.gz
```
In such cases, the samples' ID are ```E003, E008, E014```, ```<input_filename_suffix>``` will then be ```_chr22_core_K27ac_segments.bed.gz```.
For each group of samples, we need to create a file ```sample.list``` inside ```<all_cg_out_dir>/<group>``` folder. In our example, they are ```sample.list``` files in ```testdata/csrep_output/ESC``` or ```testdata/csrep_output/Brain```.

# Step 1: Adjust parameters in config file
Users will modify the ```config/config.yml``` file to adjust parameters to run CSREP. All details about what each parameter means are contained in the **Modify parameters to CSREP** section within the ```README.md``` file. To help you further understand the meaning of parameters, the section **File structure for reference in running CSREP** can help you visualize how files the input are prepared outputs are produced. 
We recommend just running snakemake using parameters that we currently sets in the config file for the example from ```testdata``` folder first. Then, the tutorial will be much easier to understand. 
 
# Step 2: Running CSREP using Snakemake 
- Please make sure that all the installation requirements are met. We recommend create a new conda environment for CSREP (refer to **Installing CSREP** section from ```README.md```)
- If you run this tutorial using your personal computer, you can run ```snakemake``` on the command line from the current working directory, which contains file ```Snakefile```.
- If you have the computing capacity of computer clusters/cloud computing, you can run ```snakemake -j 100 --cluster "qsub -V -l h_rt=4:00:00,h_data=2G"``` on the command line. Snakemake will submit jobs in parallel whenever possible, and manage your jobs. The joblog and terminal output files from each job will be saved in the home directory in your system.



