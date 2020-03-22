#----------------------------------------
#          Kaplan-Meier plot script
#----------------------------------------
#
# Supplemental R script for following study:
# Zsuzsanna Mihaly, Mate Kormos, Andras Lanczky, Magdolna Dank, Jan Budczies, A. Marcell Szusz, Balazs Gyorffy
# A meta-analysis of gene expression based biomarkers predicting outcome after tamoxifen treatment in breast cancer
# Breast Cancer Res Treat. 2013 Jul;140(2):219-32. doi: 10.1007/s10549-013-2622-y.
# 
# Input:
# - gene expression data:
#		+-----------------------------------------
#		| AffyId | Gene Expr 1 | Gene Expr 2 ...
#		+-----------------------------------------
# - clinical data
#		+-----------------------------------------
#		| AffyId | Survival time | Survival event 
#		+-----------------------------------------
#	Note: If there are other columns in the clinical table, you can specify 
#		which column has to be used for the survival data.
#	Note: The ordering of the tables has to be same.
#
# Running the script:
#	expr: expression data
#	clin: clinical data
#	event_index: column containing the survival event,
#	time_index:  column containing the survival time, 
#	affyid: if you are intrested in one Affymetrix probe ID, 
#	auto_cutoff: if this parameter is set to "true", the script finds 
#			the best cutoff value
#	quartile: if the auto_cutoff is not set, you can use the quartile option.
#			This parameter runs from 1 to 100 as a percentile, where the
#			data has to be split into the two groups. 
#			50 means the median, 25 the lower, 75 the upper quartile.
#
# Example 1:
# 	d = loadData(expr="@supplemental table 1_GEO expression data_sorted.txt", clin="@supplemental table 2_GEO clinical data_sorted.txt");
#	kmplot(d$expr, d$clin, event_index=3, time_index=4,  auto_cutoff="true")
#
# Example 2:
# 	d = loadData(expr="@supplemental table 1_GEO expression data_sorted.txt", clin="@supplemental table 2_GEO clinical data_sorted.txt");
#	kmplot(d$expr, d$clin, event_index=3, time_index=4,  affyid="213324_at", auto_cutoff="true")
#
# Example 3:
#	d = loadData(expr="expression_table.txt", clin="clinical_table.txt");
#	kmplot(d$expr, d$clin, event_index=2, time_index=3, auto_cutoff="false", quartile=50);
#
# Example 4:
# # Prepare expression data
# expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
# # Filter out low expressed genes
# # Should be more than 90% of non-zero values
# ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
# expr <- expr[, apply(expr, 2, ff)] 
# expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
# # Prepare clinical data
# clin <- mtx$merged.dat[, 1:3]
# colnames(clin)[1] <- "AffyID"
# # Run survival analysis for selected genes
# kmplot(expr, clin, event_index=2, time_index=3,  affyid = c("SND1", "MTDH"), auto_cutoff="true", transform_to_log2 = TRUE)
# # Run survival analysis for all genes
# kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)

library(survival)
library(survplot)
library(survminer)
library(cowplot)
library(gridExtra)

# demo = getParameter(c_args, "demo");
demo = "false"
if(demo == "true"){
	d = loadData();
	
	kmplot(d$expr, d$clin, auto_cutoff="true")
}


checkData = function(expr, clin){

	affyid_expr = as.character(expr[[1]]);
	affyid_clin = as.character(clin[[1]]);
	
	# check number of entries
	if(length(affyid_expr) != length(affyid_clin)){
		stop("The dimensions of data is not equal.");
	}
	
	for(i in 1:length(affyid_expr)){
		if(affyid_expr[[i]] != affyid_clin[[i]]){
			stop( paste("STOP: The ", i, "th ID of data is different.", sep="") )
		}
	}

}


auto_cutoff_surv = function(row, gene_db, time_index, event_index){

	ordered_row = order(row);
	q1 = round(length(ordered_row)*0.25);
	q3 = round(length(ordered_row)*0.75);
	
	# m = sortedrow[round(i)];
	surv = Surv(gene_db[,time_index], gene_db[,event_index]);
	
	p_values = vector(mode="numeric", length = q3-q1+1)
	min_i = 0
	min_pvalue=1
	
	for(i in q1:q3){
	
		gene_expr = vector(mode="numeric", length=length(row))
		gene_expr[ordered_row[i:length(ordered_row)]] = 1
		
		cox = summary(coxph(surv ~ gene_expr))
		
		pvalue = cox$sctest['pvalue']

		p_values[i-q1+1] = pvalue
		
		if(pvalue < min_pvalue){
			min_pvalue = pvalue
			min_i = i
		}
	
	}
	
	gene_expr = vector(mode="numeric", length=length(row))
	gene_expr[ordered_row[min_i:length(ordered_row)]] = 1
	
	
	# overwrite m (median) and gene_expr
	m = row[ordered_row[min_i]]
	
	m
}


loadData = function(exprFile="@supplemental table 1_GEO expression data_sorted.txt", clinFile="@supplemental table 2_GEO clinical data_sorted.txt"){
	
	if(file.exists(exprFile)){
		expr = read.table(exprFile, header=TRUE, sep="\t");
	}else{
		stop("The gene expression table is not exists!");
	}
	
	if(file.exists(clinFile)){
		clin = read.table(clinFile, header=TRUE, sep="\t");
	}else{
		stop("The gene expression table is not exists!");
	}
	
	list("expr"=expr, "clin"=clin);
}

mySurvplot = function(survival_data, gene_expr, xlab="Time (days)", ylab="Probability", snames = c('low', 'high'), stitle = "Expression", hr.pos=NA, use_survminer = TRUE, fileNameOut = fileNameOut) {
  # General calculations of p-value and hazard ratio
  surv <- Surv(survival_data[,1], survival_data[,2]);
  cox = summary(coxph(surv ~ gene_expr))
  pvalue = cox$sctest['pvalue'];
  hr = round(cox$conf.int[1],2)
  hr_left = round(cox$conf.int[3],2)
  hr_right = round(cox$conf.int[4],2)
  conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep=""); 
  
  # Plot survival plot
  if (use_survminer) {
    surv_expr <- data.frame(time = survival_data[, 1], status = survival_data[, 2], gene = gene_expr)
    fit <- survfit(Surv(time, status) ~ gene, data = surv_expr)
    p <- ggsurvplot(fit, data = surv_expr, risk.table = TRUE, pval = TRUE, pval.method = TRUE, conf.int = FALSE, xlab = "Time (days)", palette = c("#67A9CF", "#EF8A62"), legend = "top", legend.title = "Expression", legend.labs = c("Low", "High"))
    p <- grid.arrange(p$plot, p$table, ncol = 1, heights = c(3, 1))
    plot(p)
    # ggsave("res/text.png", plot = p$plot)
    save_plot(fileNameOut, p, base_height = 5, base_width = 5)
  } else {
    png(filename = fileNameOut)
    survplot(surv ~ gene_expr, xlab=xlab, ylab=ylab, snames = snames, stitle = stitle, hr.pos=hr.pos);
    txt = paste("HR = ", hr, conf_int, "\nlogrank P = ", signif(pvalue, 2), sep="")
    text(grconvertX(0.98, "npc"), grconvertY(.97, "npc"), labels=txt, adj=c(1, 1))
    dev.off()
  }
  # Data to return
  list(pvalue, hr, hr_left, hr_right)
}

createDirectory = function(base){
	i="";
# 	while(file.exists(paste(base, i, sep=""))){
# 		if(i==""){
# 			i=1;
# 		}else{
# 			i=i+1;
# 		}
# 	}
# }
  toDir = paste(base, i, sep="")
	
  if (!file.exists(toDir)) {
	  dir.create(toDir)
	}
  
	toDir
}


getCutoff = function(quartile, median_row, manual_cutoff="false", verbose=FALSE){

	sortedrow=order(median_row);
	minValue = median_row[sortedrow[1]]
	maxValue = median_row[sortedrow[length(sortedrow)]]

	# if manual_cutoff is true, then the user can specify a discrete Cutoff value
	# 
	if(manual_cutoff=="true"){
		m = as.numeric(quartile);
		indices = which(median_row>m)
	}else{
		quartile = as.numeric(quartile);

		if(is.na(quartile) || !is.numeric(quartile)){
			stop("The quartile parameter isn't numeric.")
		}else if(quartile<5){
			quartile = 5
		} else if(quartile > 95){
			quartile = 95
		}

		i=length(sortedrow)*quartile/100;
		m=median_row[sortedrow[round(i)]];
			
		indices = which(m<median_row)
		#sortedrow[round(i):length(sortedrow)]
	}

	if(verbose){
		print(m)
		print(minValue)
		print(maxValue)
	}

	list(m, minValue, maxValue, indices)
}

getParameter = function(c_args, id){
	res = ""
	for(i in 1:length(c_args)){

		tmp = strsplit(c_args[i], "=")
		filter_id = tmp[[1]][1]
		value = tmp[[1]][2]

		if(filter_id == paste("-", id, sep="")){
			res = value
			break
		}
	
	}
	
	res

}

kmplot = function(expr, clin, event_index=2, time_index=3, affyid="", auto_cutoff="true", quartile=50, transform_to_log2 = FALSE, cancer_type = "BRCA", fileType = "png", use_survminer = TRUE){
# checks the input: if the expression data and clinical data don't match, the script will fail.
checkData(expr, clin);

survival_data = cbind(as.numeric(clin[[time_index]]), as.numeric(clin[[event_index]]));

toDir = createDirectory("res");
resTable=rbind();

# Prepare a file for global statistics
if (!file.exists(paste0(toDir, "/global_stats.txt"))) {
  write.table( paste(c("Cancer", "Gene", "p-value", "HR", "HR_left", "HR_right", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Cutoff_type", "Cutoff_value"), collapse = "\t") , paste0(toDir, "/global_stats.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

index_arr = 2:dim(expr)[2];
if(affyid != ""){
	index = which(colnames(expr) %in% affyid);
	if(length(index) > 0){
		index_arr = c(index);
	}
}
for(j in 1:length(index_arr)){
	i=index_arr[j]
	print(paste("Processing...", colnames(expr)[i]));
	
	if (transform_to_log2) {
	  row = log2(as.numeric(expr[[i]]) + 1)
	} else {
	  row = as.numeric(expr[[i]]);
	}
	# Output expression summary statistics
	row_summary <- summary(row)
	row_summary <- data.frame(stats = names(row_summary), nums = as.vector(row_summary), stringsAsFactors = FALSE)
	print(kable(row_summary))
	
	# --------------------- CUTOFF ----------------------
	if(auto_cutoff == "true"){
	  m <- NULL
		m = auto_cutoff_surv(row, survival_data, 1, 2)
		print(paste0("Automatic cutoff at ", m))
	}else{
		# calculates lower quartile, median, or upper quartile
		tmp = getCutoff(quartile, row)
		
		m = tmp[[1]]
		minValue = tmp[[2]]
		maxValue = tmp[[3]]
		indices = tmp[[4]]
		print(paste0("Median cutoff at ", m))
	}

	# gene_expr consists 1 if m smaller then the gene expression value,
	# 0 if m bigger then the gene expression value
	gene_expr <- NULL
	gene_expr=vector(mode="numeric", length=length(row))
	gene_expr[which(row > m)] = 1 # low/high
	
	# --------------------- KMplot ----------------------
	tryCatch({
	  # draws the KM plot into a  file
	  fileNameOut <- paste0(toDir, "/", colnames(expr)[i], "_", cancer_type, ".", fileType)
	  res <- mySurvplot(survival_data, gene_expr, use_survminer = use_survminer, stitle = paste0(affyid, "\n", "Expression"), fileNameOut = fileNameOut)
		
		# Save global statistics
		pvalue = res[[1]];
		hr = res[[2]]
		hr_left = res[[3]]
		hr_right = res[[4]]
		resTable = rbind(resTable, c(pvalue, hr, hr_left, hr_right));
		write.table( paste(c(cancer_type, colnames(expr)[i], formatC(pvalue, digits = 2, format = "e"), round(c(hr, hr_left, hr_right, row_summary$nums), digits = 2), ifelse(auto_cutoff == "true", "Automatic", "Manual"), round(m, digits = 2)),  collapse = "\t") , paste0(toDir, "/global_stats.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
		
		
		
	}, interrupt = function(ex){
		cat("Interrupt during the KM draw");
		print(ex);
		dim(gene_expr);
	}, error = function(ex){
		cat("Error during the KM draw");
		print(ex);
		dim(gene_expr);
	}
	);


}

}



#' 
#' @param clin full clinical data merged with survival time and outcome. 
#' @param event_index column number to use as as outcome 
#' @param time_index column number to use for time
#' @param clinical_annotations name of the main clinical annotation category. Default: "pathologyMstage"
#' @param group1 name of the first clinical subcategory in the main catefory. Corresponds to "low" on the KM plot. Default: "m0"
#' @param group2 name of the second clinical subcategory in the main category. Corresponds to "high" on the KM plot. Default: "m1"

kmplot.clin = function(clin, event_index=2, time_index=3, clinical_annotations = "pathologyMstage", group1 = "m0", group2 = "m1", cancer_type = "BRCA", fileType = "png", use_survminer = TRUE) {
  # Full survival data
  survival_data = cbind(as.numeric(clin[[time_index]]), as.numeric(clin[[event_index]]));
  # Subset clinical annotations to the subcategories of interest
  clinical_index  <- clin[, clinical_annotations] %in% c(group1, group2) # In the main category, select all subcategories
  clinical_groups <- clin[, clinical_annotations][clinical_index] # Form a vector of labels of subcategories
  #clinical_groups <- ifelse(clinical_groups == group1, 0, 1) # Convert it to 0/1 representation of subcategories
  gene_expr <- clinical_groups # Assign to the variable traditionally used for survival analysis         
  # Subset survival data to the subcategories of interest
  survival_data <- survival_data[clinical_index, ]
  
  # Prepare output folder
  toDir = createDirectory("res");
  resTable=rbind();
  # Prepare a file for global statistics
  if (!file.exists(paste0(toDir, "/global_stats.txt"))) {
    write.table( paste(c("Cancer", "Gene", "p-value", "HR", "HR_left", "HR_right", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", paste(group1, "counts"), paste(group2, "counts")), collapse = "\t") , paste0(toDir, "/global_stats.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  # # draws the KM plot into a  file
  # if (fileType == "png") {
  #   png(paste(toDir, "/", cancer_type, "_", clinical_annotations, "_", group1, "_", group2,  ".png", sep=""));
  # }
  # if (fileType == "pdf") {
  #   pdf(paste(toDir, "/", cancer_type, "_", clinical_annotations, "_", group1, "_", group2, ".pdf", sep=""));
  # }
  # 
  # # Surv(time, event)
  # surv <- NULL
  # surv<-Surv(survival_data[,1], survival_data[,2]);
  # if (use_survminer) {
  #   res <- mySurvplot(surv, gene_expr, use_survminer = use_survminer)
  #   pvalue <- res$plot$plot_env$pval
  #   hr <- hr_left <- hr_right <- NA
  #   print(res$plot)
  # } else {
  #   res = mySurvplot(surv, gene_expr, snames = c(group1, group2), use_survminer = use_survminer)
  #   pvalue = res[[1]];
  #   hr = res[[2]]
  #   hr_left = res[[3]]
  #   hr_right = res[[4]]
  #   resTable = rbind(resTable, c(pvalue, hr, hr_left, hr_right));
  # }
  # dev.off();
  # # Save global statistics
  # write.table( paste(c(cancer_type, paste(clinical_annotations, group1, group2, sep = "-"), formatC(pvalue, digits = 2, format = "e"), round(c(hr, hr_left, hr_right), digits = 2), rep("", 6), sum(clinical_groups == group1), sum(clinical_groups == group2)),  collapse = "\t") , paste0(toDir, "/global_stats.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  
  # draws the KM plot into a  file
  fileNameOut <- paste0(toDir, "/", cancer_type, "_", clinical_annotations, "_", group1, "_", group2, ".", fileType)
  res <- mySurvplot(survival_data, gene_expr, use_survminer = use_survminer, stitle = paste0(affyid, "\n", "Expression"), fileNameOut = fileNameOut)
  
  # Save global statistics
  pvalue = res[[1]];
  hr = res[[2]]
  hr_left = res[[3]]
  hr_right = res[[4]]
  resTable = rbind(resTable, c(pvalue, hr, hr_left, hr_right));
  write.table( paste(c(cancer_type, paste(clinical_annotations, group1, group2, sep = "-"), formatC(pvalue, digits = 2, format = "e"), round(c(hr, hr_left, hr_right), digits = 2), rep("", 6), sum(clinical_groups == group1), sum(clinical_groups == group2)),  collapse = "\t") , paste0(toDir, "/global_stats.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  
  
}

