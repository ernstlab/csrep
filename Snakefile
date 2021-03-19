# Snakefile for baseline model in prediction pipeline
import pandas as pd
import os 
from itertools import product
configfile: "./config/config.yaml"
all_cg_out_dir = config["all_cg_out_dir"]
raw_segment_suffix = config['raw_segment_suffix']
raw_user_input_dir =  config['raw_user_input_dir']
training_data_folder = config["training_data_folder"]
all_ct_segment_folder = config["all_ct_segment_folder"]
sample_genome_fn = config['sample_genome_fn']
gene_reg_list = ['chr22_0']
cell_group_list = config['cell_group_list']
train_mode_list = config['train_mode_list']
sample_fn_list = list(map(lambda x: os.path.join(raw_user_input_dir, x, 'sample.list'), cell_group_list))
chrom_length_fn = config['chrom_length_fn']
num_chromHMM_state = config['chromhmm_state_num']
is_calculate_diff_two_groups = config['is_calculate_diff_two_groups'] # 0 or 1. If it is 1, which means we will calculate the differential chromatin state scores between two groups with multiple samples. If 1, the number of cell groups (cell_group_list) must be exactly 2. If 0, we will only calculate the reprentative chromatin state assignment matrix for each of the listed group in cell_group_list


def read_user_input(wildcards):
     results = []
     if is_calculate_diff_two_groups == 0: # only calculate the representative chromatin state maps for each group with cell_group_list
          fn_comb_list = list(product(cell_group_list, train_mode_list, gene_reg_list))
          results = list(map(lambda x: os.path.join(all_cg_out_dir, x[0], x[1], 'representative_data', 'average_predictions', x[2] + '_avg_pred.txt.gz'), fn_comb_list)) # x[0]: cell_group, x[1]: train_mode, x[2]: gene_reg
          return results
     elif is_calculate_diff_two_groups == 1: # we calculate the differential chromatin state map between two groups, group1 minus group2
          assert len(cell_group_list) == 2, 'Number of cell groups provided in cell_group_list should be 2 in order for csrep to calculate the differential scores between the two groups'
          output_folder = os.path.join(all_cg_out_dir, cell_group_list[0] + '_minus_' + cell_group_list[1]) # group1 minus group2 when we calcualte the differential scores
          fn_comb_list = list(product(train_mode_list, gene_reg_list))
          results = list(map(lambda x: os.path.join(output_folder, x[0], x[1] + '_avg_pred.txt.gz'), fn_comb_list)) # x[0]: train_mode, x[1]: gene_reg
          return results

rule all:
     input: 
          directory(all_ct_segment_folder), # calling rule get_all_ct_segment
          read_user_input


rule get_all_ct_segment: # this rule is called to prepare the input data for all the cell type, which will then be used to get the representative chromatin state maps for groups of sample. This rule is called only once.
     input: 
          raw_user_input_dir
     output:
          directory(all_ct_segment_folder)
     shell:
          """
          python ./scripts/get_all_ct_segment_folder.py {raw_user_input_dir} {output} {raw_segment_suffix}
          """
          
rule sample_genome: # this function will sample 10 % of genome and write those regions out into a bed file
     input:
          chrom_length_fn
     output:
          sample_genome_fn
     params:
          sample_fraction = 0.1
     shell:
          """
          python ./scripts/sample_genome_for_training.py {input} {params.sample_fraction} {output} 
          """

rule get_sample_bedfile_one_sample:
     input:
          sample_genome_fn, 
          os.path.join(raw_user_input_dir, '{ct}' + raw_segment_suffix)
     output:
          os.path.join(training_data_folder, '{ct}_train_data.bed.gz')
     params:
          output_fn_no_gz = os.path.join(training_data_folder, '{ct}_train_data.bed')
     run:
          shell("echo -e 'chrom\tstart\tend\t{wildcards.ct}' > {params.output_fn_no_gz}")
          shell("bedtools intersect -a {input[0]} -b {input[1]} -wb | awk -F'\t' 'BEGIN {{OFS=\"\t\"}} {{print $1,$2,$3,$7}}' >> {params.output_fn_no_gz}") 
          shell("gzip {params.output_fn_no_gz}")


def get_input_fn_for_average_pred_results (wildcards):
     results = []
     list_fn = os.path.join(wildcards.one_cg_out_dir, 'sample.list')
     ct_list = [line.strip() for line in open(list_fn, "r").readlines()]
     for ct in ct_list:
          results += list(map(lambda x: os.path.join(wildcards.one_cg_out_dir, 'multi_logistic', 'representative_data', 'pred_' + ct, x + '_pred_out.txt.gz'), gene_reg_list))
     return results

rule average_pred_results_csrep:
     # in order for this rule to run, the rule create_pred_multi_log_dir must run first
     input:
          get_input_fn_for_average_pred_results
     output:
          expand(os.path.join('{{one_cg_out_dir}}', "multi_logistic", 'representative_data', "average_predictions", "{gene_reg}_avg_pred.txt.gz"), gene_reg = gene_reg_list)
     params:
          list_fn = os.path.join('{one_cg_out_dir}', 'sample.list'),
          out_dir = os.path.join('{one_cg_out_dir}', 'multi_logistic', 'representative_data', 'average_predictions'),
          all_ct_pred_dir = os.path.join('{one_cg_out_dir}', 'multi_logistic', 'representative_data') # where all the subfolders of pred_<ct> are stored
     conda: 
          "config/zane_env.yml"
     shell:
          """
          python ./scripts/average_pred_results.py {params.out_dir} {params.list_fn} {all_ct_segment_folder} {params.all_ct_pred_dir} 1
          """
    
def get_training_data_one_group(wildcards):
     list_fn = os.path.join(wildcards.one_cg_out_dir, 'sample.list')
     ct_list = [line.strip() for line in open(list_fn, 'r').readlines()]
     ct_train_fn_list = list(map(lambda x: os.path.join(training_data_folder, x + '_train_data.bed.gz'), ct_list))
     return ct_train_fn_list

rule create_pred_baseline_dir:
     # no other rules need to finish first before this rule
     input: 
          get_training_data_one_group, # in order to get this, the rule get_sample_bedfile_one_sample has to be called for all the samples in this group
          all_ct_segment_folder
     params: 
          list_fn = os.path.join('{one_cg_out_dir}', 'sample.list'),
          this_predict_outDir = os.path.join('{one_cg_out_dir}', 'baseline', 'representative_data', 'average_predictions')
     output: 
          (expand(os.path.join('{{one_cg_out_dir}}', 'baseline', 'representative_data', 'average_predictions', "{gene_reg}_avg_pred.txt.gz"), gene_reg = gene_reg_list))
     conda:
          "config/zane_env.yml"
     shell:
          """
          python ./scripts/train_baseline_model.py {training_data_folder} {all_ct_segment_folder} {params.this_predict_outDir} {num_chromHMM_state} {params.list_fn}
          """
     

rule create_pred_multi_log_dir:
     input: 
          get_training_data_one_group, # in order to get this, the rule get_sample_bedfile_one_sample has to be called for all the samples in this group
          all_ct_segment_folder
     params: 
          list_fn = os.path.join('{one_cg_out_dir}', 'sample.list'),
          this_predict_outDir = os.path.join('{one_cg_out_dir}', 'multi_logistic', 'representative_data', 'pred_{train_ct}')
     output: 
          temp(expand(os.path.join('{{one_cg_out_dir}}', 'multi_logistic', 'representative_data', 'pred_{{train_ct}}', "{gene_reg}_pred_out.txt.gz"), gene_reg = gene_reg_list))
     conda:
          "config/zane_env.yml"
     shell:
          """
          python ./scripts/train_multiLog_model.py {training_data_folder} {all_ct_segment_folder} {params.this_predict_outDir} {wildcards.train_ct} {num_chromHMM_state} {params.list_fn}
          """

rule get_chrom_diff_two_group:
     input:
          expand(os.path.join(all_cg_out_dir, '{{group1}}', '{{train_mode}}', 'representative_data', 'average_predictions', '{gene_reg}_avg_pred.txt.gz'), gene_reg = gene_reg_list),
          expand(os.path.join(all_cg_out_dir, '{{group2}}', '{{train_mode}}', 'representative_data', 'average_predictions', '{gene_reg}_avg_pred.txt.gz'), gene_reg = gene_reg_list),
     output:
          expand(os.path.join(all_cg_out_dir, '{{group1}}_minus_{{group2}}', '{{train_mode}}', '{gene_reg}_avg_pred.txt.gz'), gene_reg = gene_reg_list)
     params:
          group1_dir = os.path.join(all_cg_out_dir, '{group1}', '{train_mode}', 'representative_data', 'average_predictions'),
          group2_dir = os.path.join(all_cg_out_dir, '{group2}', '{train_mode}', 'representative_data', 'average_predictions'),
          out_dir = os.path.join(all_cg_out_dir, '{group1}_minus_{group2}', '{train_mode}'),
     shell:
          """
          python ./scripts/calculate_diff_rep_state_map.py {params.group1_dir} {params.group2_dir} {params.out_dir} {num_chromHMM_state}
          """
