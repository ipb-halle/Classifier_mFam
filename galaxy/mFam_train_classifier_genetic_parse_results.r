#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)



# ---------- Parse results and print statistics ----------
gen_mFam_results_files <- list.files(".", pattern="*", recursive=TRUE, full.names=TRUE)
gen_mFam_results_files <- gen_mFam_results_files[grep(".*Results\\.tsv", gen_mFam_results_files, invert=FALSE)]

gen_mFam_result_table <- data.frame(filename=0, score=0)
for (i in gen_mFam_results_files) {
	gen_mFam_results <- read.table(file=i, sep='\t', header=TRUE, stringsAsFactors=TRUE, fill=TRUE)
	
	score_auc_pr <- sum(gen_mFam_results$AUC.PR, na.rm=TRUE)
	score_num_classes <- nrow(gen_mFam_results)
	score_tpr <- sum(gen_mFam_results$TPR.for.FPR...5., na.rm=TRUE)
	
	# Calculate fitness score
	if (! is.numeric(score_auc_pr)) score_auc_pr <- 0
	if (! is.numeric(score_tpr)) score_tpr <- 0
	if ((! is.numeric(score_num_classes)) | (score_num_classes < 1)) score_num_classes <- 1
	fitness_score <- score_auc_pr / score_num_classes
	
	print(paste0("GA file ", i, ", Fitness Score: ", fitness_score))
	gen_mFam_result_table[nrow(gen_mFam_result_table)+1, "filename"] <- as.character(i)
	gen_mFam_result_table[nrow(gen_mFam_result_table), "score"] <- as.numeric(fitness_score)
}

gen_mFam_results_max <- read.table(file=gen_mFam_result_table$filename[which.max(gen_mFam_result_table$score)], sep='\t', header=TRUE, stringsAsFactors=TRUE, fill=TRUE)
print("")
print(paste0("Compound library with maximum score: ", gen_mFam_result_table$filename[which.max(gen_mFam_result_table$score)]))
print(gen_mFam_results_max[, c("Substance.class", "AUC.PR")])
