log_folder=$(readlink -f $1)

for f in ${log_folder}/snakejob.*.e*
do
	num_timestamp_line=$(cat $f | grep '^\Finished' | wc -l)
	if [ $num_timestamp_line != 1 ]
	then
		job=$(echo $f | awk -F'/' '{print $NF}' | awk -F'.' '{print $2}')
		out_f=$(echo $f | tr '.e' '.o')
		echo Error file associated with an unfinished/failed job: $f
		echo Check this file and the accompanying terminal output file for error message: $out_f
		echo Asosciated job: $job
		echo
	fi
done