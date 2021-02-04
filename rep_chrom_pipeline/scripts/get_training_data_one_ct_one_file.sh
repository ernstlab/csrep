##### GET COMMAND LINE ARGUMENTS #####
cell_group=$1 # The only command line argument is the cell type that we are looking at to do more analysis
######################################
num_chromHMM_state=18
sample_genome_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/sample_genome.py
roadmap_outDir=/u/home/h/havu73/project-ernst/diff_pete/roadmap
posterior_all_ct_dir=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior
training_data_outDir=${roadmap_outDir}/training_data/
mkdir -p $training_data_outDir # create directory for this cell cell group
roadmap18_segmentation_folder=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/segmentation/hg19/
ct_list_fn=${roadmap_outDir}/${cell_group}.list

# set up code file
shDir=/u/home/h/havu73/project-ernst//source_pete/jobs
shPrefix=get_sample_data_${cell_group}
sh_fn=${shDir}/${shPrefix}.sh
rm -f $sh_fn

# first, we have to sample regions on the genome that we will use to train the model
sample_bin_fn=${roadmap_outDir}/sample_genome_regions.gz
command="python $sample_genome_code $roadmap18_segmentation_folder $ct_list_fn $sample_bin_fn"
# $command

# second, we will use bedtools to filter the segmentation of sampled regions in all cell types of interest
ct_list=""
while read line;  
do
        ct_list="$ct_list $line"
done < $ct_list_fn
for ct in $ct_list
do
        ct_segment_fn=$roadmap18_segmentation_folder/$ct/${ct}_18_core_K27ac_segments.bed.gz
        ct_sample_segment_fn=${training_data_outDir}/${ct}_train_data.bed
        if [ ! -f ${ct_sample_segment_fn}.gz ]
        then 
			echo "rm -f ${ct_sample_segment_fn}.gz">> $sh_fn
			echo "echo -e 'chrom\tstart\tend\t${ct}' >${ct_sample_segment_fn}">> $sh_fn
			echo "bedtools intersect -a $sample_bin_fn -b $ct_segment_fn -wb | awk -F'\t' 'BEGIN {OFS=\"\t\"} {print \$1,\$2,\$3,\$7}' >> $ct_sample_segment_fn" >> $sh_fn 
			echo "gzip ${ct_sample_segment_fn}" >> $sh_fn
			echo Done with cell type: $ct
	    fi
done



