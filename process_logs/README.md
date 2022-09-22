# Why do I need this file
This folder contains code that you can use to check for unfinished jobs' log files in the event that your snakemake run ended without finishing all the jobs. Usually, when you use snakemake to submit jobs in the CSREP pipeline onto your computing cluster, snakemake will handle all "job-sitting" work for you. Sometimes, due to file system latency (a job is finished but all the output files are not reported as being completely written out, due to system latency), snakemake will crash. In this case, you can wait for a few minutes or until all submitted jobs on your computing cloud are finished, and run the script ```find_unfinished_snakemake_jobs.sh```  to find a list of jobs that are unfinished/failed.

# What to do when snakemake crashes
- First, wait for a few minutes or until all jobs submitted by snakemake are finished. 
- Second, ```snakemake -n --snakefile <path/to/Snakefile>``` to see all the jobs that are about to be rerun, should we rerun the pipeline again. Snakemake will know to pick up from unfinished jobs and not repeat finished jobs. 
- Third, you can consider running snakemake again to submit unfinished jobs, which can resolve problems associated with file system latency. 
- Fourth, you can also run ```find_unfinished_snakemake_jobs.sh``` to find a list of unfinished jobs and their log files, which contains error messages of the respective jobs. 

# How do I run ```find_unfinished_snakemake_jobs.sh```
- If you run snakemake directly on a computer, then you will not need code from this directory. Because running snakemake on a single machine means running jobs sequentially, and error messages will appear should there are any, and snakemake will crash. 
- If you run snakemake using the computing cluster where jobs corresponding to different components of the pipeline can be submitted in parallel, you can use a command such as ```snakemake -j 100 --cluster "qsub -V -l h_rt=4:00:00,h_data=2G -o <path/to/log/folder> -e <path/to/log/folder>" --snakefile <path/to/snakefile>```. Then, ```<path/to/log/folder>``` show the local path to a folder where you save all the log files returned by individual jobs on the computing cluster. **Then, you can use the command ```find_unfinished_snakemake_jobs.sh <path/to/log/folder>``` to search within this folder for log files corresponding to failed jobs. **
- If there are no terminal output from running the above command, that means all jobs finished successfully. 
- Then, you can look into these log files for the error messages. 
- If you need assistance, please reach out to us (jason dot ernst at ucla dot edu and havu73 at ucla dot edu). 
- Note: between each run of snakemake, you can remove all log files from the folder ```<path/to/log/folder>```, to keep track of those associated with current run of snakemake.

