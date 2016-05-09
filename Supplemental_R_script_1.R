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

library(survival)
library(survplot)


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


auto_cutoff = function(row, gene_db, time_index, event_index){

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

mySurvplot = function(surv, gene_expr, xlab="Time (years)", ylab="Probability", snames = c('low', 'high'), stitle = "Expression", hr.pos=NA){
	survplot(surv ~ gene_expr, xlab=xlab, ylab=ylab, snames = snames, stitle = stitle, hr.pos=hr.pos);
	
	cox = summary(coxph(surv ~ gene_expr))
	
	pvalue=cox$sctest['pvalue'];
	hr = round(cox$conf.int[1],2)
	hr_left = round(cox$conf.int[3],2)
	hr_right = round(cox$conf.int[4],2)
	
	conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep=""); 
	
	txt = paste("HR = ", hr, conf_int, "\nlogrank P = ", signif(pvalue, 2), sep="")
	text(grconvertX(0.98, "npc"), grconvertY(.97, "npc"),
         labels=txt,
         adj=c(1, 1))
	
	list(pvalue, hr, hr_left, hr_right)
}

createDirectory = function(base){
	i="";
	while(file.exists(paste(base, i, sep=""))){
		if(i==""){
			i=1;
		}else{
			i=i+1;
		}
	}
	toDir = paste(base, i, sep="")
	dir.create(toDir)
	
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

kmplot = function(expr, clin, event_index=3, time_index=4, affyid="", auto_cutoff="true", quartile=50){

# checks the input: if the expression data and clinical data don't match, the script will fail.
checkData(expr, clin);

survival_data = cbind(as.numeric(clin[[time_index]]), as.numeric(clin[[event_index]]));

toDir = createDirectory("res");
resTable=rbind();

index_arr = 2:dim(expr)[2];
if(affyid != ""){
	index = which(colnames(expr) == affyid);
	if(length(index) > 0){
		index_arr = c(index);
	}
}
for(j in 1:length(index_arr)){
	i=index_arr[j]
	print(paste("Processing...", colnames(expr)[i]));
	
	row = as.numeric(expr[[i]]);
	# Output expression summary statistics
	row_summary <- summary(row)
	data.frame(stats = names(row_summary), nums = as.vector(row_summary), stringsAsFactors = FALSE) %>% kable %>% print
	
	# --------------------- CUTOFF ----------------------
	if(auto_cutoff == "true"){
		m = auto_cutoff(row, survival_data, 1, 2)
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
	gene_expr=vector(mode="numeric", length=length(row))
	gene_expr[which(m < row)] = 1
	
	# --------------------- KMplot ----------------------
	tryCatch({
		# draws the KM plot into a png file
		
		png(paste(toDir, "/", colnames(expr)[i], ".png", sep=""));

		# Surv(time, event)
		surv<-Surv(survival_data[,1], survival_data[,2]);

		res = mySurvplot(surv, gene_expr)
		pvalue = res[[1]];
		hr = res[[2]]
		hr_left = res[[3]]
		hr_right = res[[4]]
		
		resTable = rbind(resTable, c(pvalue, hr, hr_left, hr_right));

		dev.off();
		
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




