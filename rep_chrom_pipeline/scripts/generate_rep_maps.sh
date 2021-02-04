	##### GET COMMANDLINE ARGUMENTS #####
cell_group=$1
train_mode=$2 # ['baseline', 'multi_logistic', 'chromHMM_posterior']
#####################################
##### CODE FILES ####
pred_baseline_code=/u/home/h/havu73/project-ernst/source_pete/rep_chrom_pipeline/scripts/train_baseline_model.py
pred_multiLog_code=/u/home/h/havu73/project-ernst/source_pete/rep_chrom_pipeline/scripts/train_multiLog_model.py
calculate_avg_pred_code=/u/home/h/havu73/project-ernst/source_pete/rep_chrom_pipeline/scripts/average_pred_results.py
#####################
state_annot_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/state_annotation.txt # annotation of the 18 states in ROADMAP 18-state models
roadmap_outDir=/u/home/h/havu73/project-ernst/diff_pete/roadmap/
cg_outDir=$roadmap_outDir/${cell_group}/${train_mode}
mkdir -p $cg_outDir
cell_group_fn=${roadmap_outDir}/${cell_group}.list
all_ct_segment_folder=/u/home/h/havu73/project-ernst/diff_pete/roadmap/all_ct_segments
training_data_folder=/u/home/h/havu73/project-ernst/diff_pete/roadmap/training_data
num_chromHMM_state=18
shDir=/u/home/h/havu73/project-ernst/source_pete/jobs
#rm -f $draw_job_fn
job_index=1
# first, we will get all the cell types related to ${cell_group}
ct_list=""
while read line;  
do
	ct_list="$ct_list $line"
done < $cell_group_fn

# now create the validation jobs for each cell type as the validation cell type
for ct in $ct_list
do
	sh_file_name=${shDir}/pred_${cell_group}_${job_index}.sh
	rm -f $sh_file_name # so that we can create new code file
	this_ct_outdir=${cg_outDir}/pred_${ct}/
	mkdir -p $this_ct_outdir
	num_file_in_this_folder=$(ls $this_ct_outdir/*.gz | wc -l | awk '{print $1}')
	if [[ ${num_file_in_this_folder} != 316 ]] 
	then
		echo "cell_group: $cell_group . Cell_type: $ct . validation of this cell type is not available"
		if [ $train_mode == "baseline" ] 
		then
			command="python $pred_baseline_code ${training_data_folder} ${all_ct_segment_folder} ${this_ct_outdir} $ct ${num_chromHMM_state} $cell_group_fn"
			echo $command > $sh_file_name
		fi
		if [ $train_mode == 'multi_logistic' ]
		then
			command="python $pred_multiLog_code ${training_data_folder} ${all_ct_segment_folder} ${this_ct_outdir} $ct ${num_chromHMM_state} $cell_group_fn"
			echo $command >$sh_file_name
		fi
		chmod +x $sh_file_name
		job_index=$(($job_index + 1))
		echo Done with cell type: $ct
	fi
done

this_cg_avg_outdir=${cg_outDir}/average_predictions/
mkdir -p ${this_cg_avg_outdir}
num_file_in_this_folder=$(ls $this_cg_avg_outdir | wc -l | awk '{print $1}')
if [[ ${num_file_in_this_folder} != 316 ]]
then
	command="python $calculate_avg_pred_code $this_cg_avg_outdir $cell_group_fn $all_ct_segment_folder ${cg_outDir}"
	echo $command >${shDir}/average_${cell_group}.sh 
	chmod +x ${shDir}/average_${cell_group}.sh 
fi
