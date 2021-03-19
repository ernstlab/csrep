# Snakefile for baseline model in prediction pipeline
import pandas as pd
import os 
configfile: "./config/config.yaml"
roadmap18_segmentation_folder='../../data/roadmap_epigenome/18_core_K27ac_model_downloaded/segmentation/hg19/'
all_cg_out_dir = config["all_cg_out_dir"]
all_ct_chrom_pos_folder = config['all_ct_chrom_pos_folder']
gene_reg_fn = config["gene_region_fn"]
gene_reg_list = ['chrX_9', 'chr1_0']#[line.strip() for line in open(gene_reg_fn, "r").readlines()]
train_mode_list = ['baseline']#['baseline', 'multi_logistic', 'chromHMM_posterior']
cell_group_list = ['iPSC', 'Heart', 'Muscle', 'Skin', 'Sm_Muscle', 'blood', 'Blood_other', 'ESC', 'Brain', 'ES-deriv'] 
#train_mode = config["train_mode"]
training_data_folder = config["training_data_folder"]
all_ct_segment_folder = config["all_ct_segment_folder"]
state18_annot_fn = config['state18_annot_fn']
num_chromHMM_state = config['chromhmm_state_num']
num_score_bins = config['num_score_bins_roc']


rule all:
     input: 
          expand(os.path.join(all_cg_out_dir, '{cell_group}', "{train_mode}", 'representative_data', "average_predictions", "{gene_reg}_avg_pred.txt.gz"), cell_group = cell_group_list, gene_reg = gene_reg_list, train_mode = train_mode_list),


def get_input_fn_for_average_pred_results (wildcards):
     results = []
     list_fn = os.path.join(wildcards.one_cg_out_dir, 'sample.list')
     ct_list = [line.strip() for line in open(list_fn, "r").readlines()]
     for ct in ct_list:
          results += list(map(lambda x: os.path.join(wildcards.one_cg_out_dir, 'multi_logistic', 'representative_data', 'pred_' + ct, x + '_pred_out.txt.gz'), gene_reg_list))
     return results

rule sample_genome: 
     input:
          roadmap18_segmentation_folder,
          os.path.join(all_cg_out_dir, '{cg}', 'sample.list')
     output:
          os.path.join(training_data_folder)
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
          
rule create_pred_baseline_dir:
     # no other rules need to finish first before this rule
     input: 
          training_data_folder,
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
          training_data_folder,
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

rule create_pred_chromHMM_pos_dir:
     input:
          all_ct_chrom_pos_folder, 
     params:
          this_predict_outDir = os.path.join('{one_cg_out_dir}', 'chromHMM_posterior', 'representative_data', 'average_predictions'),
          list_fn = os.path.join('{one_cg_out_dir}', 'sample.list'),
     output: 
          (expand(os.path.join('{{one_cg_out_dir}}', 'chromHMM_posterior', 'representative_data', 'average_predictions', "{gene_reg}_avg_pred.txt.gz"), gene_reg = gene_reg_list))
     conda:
          "config/zane_env.yml"
     shell:
          """
          python ./scripts/train_chromHMM_posterior_model.py {all_ct_chrom_pos_folder} {params.this_predict_outDir} {num_chromHMM_state} {params.list_fn}
          """
