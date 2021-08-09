library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library('ggpubr')
read_confusion_matrix_df <- function(fn, train_mode){
  df <- read.table(fn, sep = '\t', header = T, stringsAsFactors = F)
  colnames(df) <- c('ct', colnames(df)[-1])
  df$tm <- rep(train_mode, nrow(df))
  return(df)
}

get_avg_confusion_matrix_fn <- function(train_mode, cell_group_dir){
	return(file.path(cell_group_dir, train_mode, 'average_roc_all_ct/summary_confusion_matrix.txt.gz'))
}

draw_roc_plot_compare_train_mode <- function(plot_df, out_dir, state, cell_group){
	save_fn <- file.path(out_dir, paste0("avg_roc_state", state, ".png"))
	plot_title <- paste0(cell_group, "_diag_conf_S", state)
	p <- ggplot(data = plot_df, aes_string(x = 'tm', y = paste0('fraction_true_S', state,'_correct'), group = 'ct')) +
	geom_line(color = 'red') +
	ylim(0,1) +
	theme_bw() +
	ggtitle(plot_title)+
	theme(plot.title = element_text(size = 6, face = 'bold'),
	      axis.title = element_blank(),
	      axis.text.x = element_text(size = 4.5),
	      axis.text.y = element_text(size = 4.5),
	      legend.position = 'None')
	return(p)
}

arrange_18_state_plots <- function(plot_list, out_dir){
	# I could not figure out how to make this function work for variable number of states, so here we fix it to 18 states. If the number of states are different from 18, then this code will crash
	ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]],  ncol = 3, nrow = 6)
	save_fn <- file.path(out_dir, 'diag_confusion_matrix_combined_states.png')
	ggsave(save_fn, width = 7.5, height = 15)
}

main <- function(num_train_mode, cell_group_dir, out_dir, train_mode_list, num_chromHMM_state, cell_group){
	# 1. get the list of files that contain the data of average roc for different training modes
	roc_fn_list <- sapply(train_mode_list, get_avg_confusion_matrix_fn, cell_group_dir = cell_group_dir) 
	roc_fn_list <- as.character(roc_fn_list)
	# 2. get the list of dataframs from the list of file names that we just got
	roc_df_list <- list()
	for (tm_index in seq(num_train_mode)){
		this_tm_df <- read_confusion_matrix_df(roc_fn_list[tm_index], train_mode_list[tm_index])
		roc_df_list[[tm_index]] <- this_tm_df
	}
	# loop through each state
	plot_list <- list()
	for (state in 1:num_chromHMM_state){
		columns_to_select <- c('ct', paste0('fraction_true_S', state,'_correct'), 'tm')
		state_confusion_diag_df_list <- lapply(roc_df_list, FUN = function(x) (x %>% select(columns_to_select))) # for each df, select only the ones associated with the state we are trying to create a plot for
		plot_df <- bind_rows(state_confusion_diag_df_list)
		print(plot_df)
		print("before draw_roc_plot_compare_train_mode")
		p <- draw_roc_plot_compare_train_mode(plot_df, out_dir, state, cell_group)
		plot_list[[state]] <- p 
		print(paste0("Done with drawing plots for state", state))
	}
	print('Going to call a function that only works when num_chromHMM_state = 18. Otherwise, it will crash!')
	arrange_18_state_plots(plot_list, out_dir)
	print("Done combining all the states' plots")
}


args = commandArgs(trailingOnly=TRUE)
NUM_SURE_ARGS = 5
if (length(args) < NUM_SURE_ARGS)
{
	stop("wrong command line argument format", call.=FALSE)
}
num_train_mode <- as.integer(args[1]) # output of ChromHMM OverlapEnrichment, also input to this program
if(length(args) != (NUM_SURE_ARGS + num_train_mode)){
	stop("Number of arguments do not match the number of specified train_mode", call.=FALSE)
}
cell_group_dir <- args[2] # where the figures should be stored
if (! dir.exists(cell_group_dir)){
	stop(paste0("cell_group_dir DOES NOT EXIST: ", cell_group_dir), call.=FALSE)
}
out_dir <- args[3] # name of the enrichment context, so that it can be title of the figure
dir.create(out_dir, recursive = TRUE)
num_chromHMM_state <- as.integer(args[4])
cell_group <- args[5]
train_mode_list <- args[(NUM_SURE_ARGS + 1):length(args)]
main(num_train_mode, cell_group_dir, out_dir, train_mode_list, num_chromHMM_state, cell_group)