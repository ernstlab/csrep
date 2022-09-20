This folder contains some useful scripts for your to preprocess or post-process data. You can refer to the docstrings provided in each python script if you need further information about individual scripts. 
```get_chrom_length_from_segmentation.sh <input_segment_fn> <output_fn> <chrom_list_fn>```: Given an input sample's chromatin state segmentation file, this script will find the lengths of all main chromosomes. The output will include 3 columns corresponding to ```chromosomes	0	length```. Note that the output file follows a bed format. The third argument ```chrom_list_fn``` should specify the chromosomes for which you would like to get the length, with each line showing one chromosome (1, 2, etc., X). Among these chromosomes, those without segmentation data will output length 0 (example: some input segmentation may not contain any data for chromosome Y). We already provide the files ```chrom.human``` and ```chrom.mouse``` that specify the lists of chromosomes in human and mouse genome. The use case for this script, i.e. reason why we want to extract the chromsome length from segmentation data is so that we can generate a list of randomly-selected regions that consitute the training data of CSREP. The output of this script may be used for ```../scripts/sample_genome_for_training.py``` 
An example call of this code is ```./get_chrom_length_from_segmentation.sh ../testdata/raw_data/E003/E003_chr22_core_K27ac_segments.bed.gz ./trial.output ./chrom.human ```. This will produce the ```trial.output``` as follows: 

```
chr1    0       
chr2    0       
chr3    0       
chr4    0       
chr5    0       
chr6    0       
chr7    0       
chr8    0       
chr9    0       
chr10   0       
chr11   0       
chr12   0       
chr13   0       
chr14   0       
chr15   0       
chr16   0       
chr17   0       
chr18   0       
chr19   0       
chr20   0       
chr21   0       
chr22   0       51304400
chrX    0       
chrY    0       
```
The reason for this is because currently, the provided data inside ```../testdata/raw_data``` only contains segmentation data in chr22, therefore, for all other chromosomes, the output length is left blank. The code ```../scripts/sample_genome_for_training.py``` will be able to read this file and generate the appropriate sample genomic bins. 

```python compress_segmentation_data.py --segment_bed_fn <segment_bed_fn> --output_fn <output_fn>```: This code takes in a segmentation file and reaarranges the rows such that if consecutive genomic bins are annotated as the same state, they will be combined into one line. At the same time, the output file will be sorted by chrom, start_bp such that it will be readily available to use for downstream analysis with bedtools or pybedtools. 
```helper.py```: This file contains some useful functions for file-managements, reading input data, etc. that can be used by other scripts. 