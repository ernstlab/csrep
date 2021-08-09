library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library('pheatmap')

get_raw_state_number_from_conf_names <- function(state_name){
	# change from 'true_E1'. to 1
 	state_name <- (unlist(strsplit(state_name, "_"))[2]) # E1
 	state_name <- substr(state_name, 2, nchar(state_name)) # ex: '1' --> '100'
 	return (as.numeric(as.vector(state_name)))
}

get_confusion_df_with_annotated_state_names <- function(state_annot_fn, confusion_fn){
	confusion_df <- read.table(confusion_fn, header = TRUE, sep = '\t')
	state_annot_df <- read.table(state_annot_fn, header = TRUE, sep = '\t') # df of state anontations such as colors and mnenomics
	confusion_df$X <- as.character(confusion_df$X)
	conf_state_df <- data.frame(state = sapply(confusion_df$X, get_raw_state_number_from_conf_names)) # we need as character because soemhow R can inteprete a vector of 'true_E1', 'true_E2', etc. as integer. R is so stupid
	state_annot_df <- left_join(conf_state_df, state_annot_df) # join on state
	state_annot_df['long_state_name'] <- paste0(state_annot_df$state, '_', state_annot_df$mnenomic)
	rownames(confusion_df) <- state_annot_df$long_state_name
	confusion_df <- confusion_df %>% select(-'X') 
	return(confusion_df)
}

draw_confusion_matrix <- function(confusion_fn, save_fn, state_annot_fn){
	confusion_df <- get_confusion_df_with_annotated_state_names(state_annot_fn, confusion_fn)
	pheatmap(confusion_df, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = TRUE, show_rownames = TRUE, angle_col = 90, filename = save_fn)
	return(confusion_df)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3)
{
	stop("wrong command line argument format", call.=FALSE)
}
args = commandArgs(trailingOnly=TRUE)
confusion_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
save_fn <- args[2] # where the figures should be stored
state_annot_fn <- args[3] # where data of states' mnenomics and colors are stored
confusion_confusion_df <- draw_confusion_matrix(confusion_fn, save_fn, state_annot_fn) 
