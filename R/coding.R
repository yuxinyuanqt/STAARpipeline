coding <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                   rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                   SPA_p_filter=FALSE,p_filter_cutoff=0.05,use_ancestry_informed=FALSE,find_weight=FALSE,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)
	n_pheno <- obj_nullmodel$n.pheno

	## SPA status
	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}

	## get SNV id, position, REF, ALT (whole genome)
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	position <- as.numeric(seqGetData(genofile, "position"))
	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	rm(position)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	################################################
	#           Coding
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.coding]
	
	if(length(variant.id.gene)==0)
	{
	  results_coding <- list(plof = c(),
	                         plof_ds = c(),
	                         missense = c(),
	                         disruptive_missense = c(),
	                         synonymous = c())
	  seqResetFilter(genofile)
	  return(results_coding)
	}

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	## Annotation
	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

	if(variant_type=="SNV")
	{
		if(Use_annotation_weights)
		{
			for(k in 1:length(Annotation_name))
			{
				if(Annotation_name[k]%in%Annotation_name_catalog$name)
				{
					Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
					Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

					if(Annotation_name[k]=="CADD")
					{
						Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
					}

					if(Annotation_name[k]=="aPC.LocalDiversity")
					{
						Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
						Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
						Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
					}
					Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
				}
			}

			Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
			colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
		}
	}
	
	## get AF, Missing rate
	AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
	REF_AF <- AF_AC_Missing$af
	Missing_rate <- AF_AC_Missing$miss
	variant_maf_cutoff_filter <- rare_maf_cutoff

	################################################
	#                  plof_ds
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	results_information <- Genotype_sp$results_information
	rm(Genotype_sp)
	gc()
	
	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  MAF.in <- results_information$MAF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  rm(results_information)
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }
	  
	  if(use_ancestry_informed == TRUE)
	  {
	    Geno <- as.matrix(Geno)
	  }
	} else
	{
	  Geno <- NULL
	  MAF.in <- numeric()
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
				pvalues_plof_ds <- pvalues
			}
		}else{
			try(pvalues <- STAAR_Binary_SPA_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
		}
	}else
	{
		try(pvalues <- MultiSTAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	results_plof_ds <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "plof_ds"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results_plof_ds <- rbind(results_plof_ds,results_temp)
	}

	if(!is.null(results_plof_ds))
	{
		if(!use_SPA)
		{
			colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
			colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O","STAAR-O")
		}else
		{
			colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
			colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_plof_ds)[dim(results_plof_ds)[2]] <- c("STAAR-B")

		}
	}

	#####################################################
	#                      plof
	#####################################################
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	results_information <- Genotype_sp$results_information
	rm(Genotype_sp)
	gc()
	
	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  MAF.in <- results_information$MAF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  rm(results_information)
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }
	  
	  if(use_ancestry_informed == TRUE)
	  {
	    Geno <- as.matrix(Geno)
	  }
	} else
	{
	  Geno <- NULL
	  MAF.in <- numeric()
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
				pvalues_plof <- pvalues
			}
		}else{
			try(pvalues <- STAAR_Binary_SPA_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
		}
	}else
	{
		try(pvalues <- MultiSTAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	results_plof <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "plof"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results_plof <- rbind(results_plof,results_temp)
	}

	if(!is.null(results_plof))
	{
		if(!use_SPA)
		{
			colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
			colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O","STAAR-O")
		}else
		{
			colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
			colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_plof)[dim(results_plof)[2]] <- c("STAAR-B")
		}
	}

	#############################################
	#             synonymous
	#############################################
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.synonymous],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.synonymous],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	results_information <- Genotype_sp$results_information
	rm(Genotype_sp)
	gc()
	
	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  MAF.in <- results_information$MAF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  rm(results_information)
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }
	  
	  if(use_ancestry_informed == TRUE)
	  {
	    Geno <- as.matrix(Geno)
	  }
	} else
	{
	  Geno <- NULL
	  MAF.in <- numeric()
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
				pvalues_synonymous <- pvalues
			}
		}else{
			try(pvalues <- STAAR_Binary_SPA_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
		}
	}else
	{
		try(pvalues <- MultiSTAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	results_synonymous <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "synonymous"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results_synonymous <- rbind(results_synonymous,results_temp)
	}

	if(!is.null(results_synonymous))
	{
		if(!use_SPA)
		{
			colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
			colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O","STAAR-O")
		}else
		{
			colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
			colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results_synonymous)[dim(results_synonymous)[2]] <- c("STAAR-B")
		}

	}

	#################################################
	#        missense
	#################################################
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.missense]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.missense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.missense],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	results_information <- Genotype_sp$results_information
	rm(Genotype_sp)
	gc()
	
	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  MAF.in <- results_information$MAF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  rm(results_information)
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }
	  
	  if(use_ancestry_informed == TRUE)
	  {
	    Geno <- as.matrix(Geno)
	  }
	} else
	{
	  Geno <- NULL
	  MAF.in <- numeric()
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
				pvalues_m <- pvalues
			}
		}else{
			try(pvalues <- STAAR_Binary_SPA_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
		}
	}else
	{
		try(pvalues <- MultiSTAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	results <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results <- rbind(results,results_temp)
	}

	#################################################
	#         disruptive missense
	#################################################
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.dmissense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.dmissense],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	results_information <- Genotype_sp$results_information
	rm(Genotype_sp)
	gc()
	
	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  MAF.in <- results_information$MAF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  rm(results_information)
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }
	  
	  if(use_ancestry_informed == TRUE)
	  {
	    Geno <- as.matrix(Geno)
	  }
	} else
	{
	  Geno <- NULL
	  MAF.in <- numeric()
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
				pvalues_ds <- pvalues
			}
		}else{
			try(pvalues <- STAAR_Binary_SPA_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
		}
	}else
	{
		try(pvalues <- MultiSTAAR_sp(Geno,MAF.in,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "disruptive_missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results <- rbind(results,results_temp)
	}

	if(!is.null(results))
	{
		if(!use_SPA)
		{
			colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
			colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
		}else
		{
			colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
			colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results)[dim(results)[2]] <- c("STAAR-B")
		}

		if(dim(results)[1]==1)
		{
			if(results[3]!="disruptive_missense")
			{
				if(!use_SPA)
				{
					results <- cbind(results,matrix(1,1,6))
					colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
					results_missense <- results
					results_ds <- c()
				}else{
					results <- cbind(results,matrix(1,1,2))
					colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
					results_missense <- results
					results_ds <- c()
				}

				if(use_ancestry_informed == TRUE & find_weight == TRUE & !use_SPA){
					results_weight <- results_weight1 <- results_weight2 <- c()
					for(i in 1:ncol(pvalues_m$results_weight)){
						results_m_weight <- pvalues_m$results_weight[-c(1,2),i]
						results_m_weight <- unlist(pvalues_m$results_weight[,i][c(5:length(pvalues_m$results_weight[,i]), 4,3)])
						names(results_m_weight) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]

						results_m_weight <- c(results_m_weight, rep(1,6))
						names(results_m_weight)[(length(results_m_weight)-5):length(results_m_weight)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight <- cbind(results_weight, results_m_weight)
						colnames(results_weight)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues_m$results_weight1)){
						results_m_weight1 <- pvalues_m$results_weight1[-c(1,2),i]
						results_m_weight1 <- unlist(pvalues_m$results_weight1[,i][c(5:length(pvalues_m$results_weight1[,i]), 4,3)])
						names(results_m_weight1) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]

						results_m_weight1 <- c(results_m_weight1, rep(1,6))
						names(results_m_weight1)[(length(results_m_weight1)-5):length(results_m_weight1)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight1 <- cbind(results_weight1, results_m_weight1)
						colnames(results_weight1)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues_m$results_weight2)){
						results_m_weight2 <- pvalues_m$results_weight2[-c(1,2),i]
						results_m_weight2 <- unlist(pvalues_m$results_weight2[,i][c(5:length(pvalues_m$results_weight2[,i]), 4,3)])
						names(results_m_weight2) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]

						results_m_weight2 <- c(results_m_weight2, rep(1,6))
						names(results_m_weight2)[(length(results_m_weight2)-5):length(results_m_weight2)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight2 <- cbind(results_weight2, results_m_weight2)
						colnames(results_weight2)[i] <- c(i-1)
					}

					results_weight_m <- results_weight
					results_weight1_m <- results_weight1
					results_weight2_m <- results_weight2

					rownames(pvalues_m$weight_all_1) <-rownames(pvalues_m$weight_all_2) <- unique(obj_nullmodel$pop.groups)

					results <- list(results,
					                weight_all_1 = pvalues_m$weight_all_1,
					                weight_all_2 = pvalues_m$weight_all_2,
					                results_weight = results_weight_m,
					                results_weight1 = results_weight1_m,
					                results_weight2 = results_weight2_m)
				}
			}else
			{
				results_missense <- c()
				results_ds <- results
				results <- c()
			}
		}
		
		if(!is.null(results))
		{
			if(dim(results)[1]==2)
			{
				if(!use_SPA)
				{
					results_m <- c(results[1,],rep(0,6))
					names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
					results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT(1,25)","SKAT(1,1)","Burden(1,25)","Burden(1,1)","ACAT-V(1,25)","ACAT-V(1,1)")]
					apc_num <- (length(results_m)-19)/6
					p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
					results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
					results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num,6*apc_num+9)]))
					results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
					results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
					results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
					results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
					results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

					results_ds <- c()
					results_ds <- rbind(results_ds,results[2,])

					results <- c()
					results <- rbind(results,results_m)

					#results_weight
					if(use_ancestry_informed == TRUE & find_weight == TRUE){
						results_weight_m <- results_weight1_m <- results_weight2_m <- c()

						for(i in 1:ncol(pvalues_m$results_weight)){
							results_m_weight <- pvalues_m$results_weight[-c(1,2),i]
							results_m_weight <- unlist(pvalues_m$results_weight[,i][c(5:length(pvalues_m$results_weight[,i]), 4,3)])
							names(results_m_weight) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

							results_m_weight <- c(results_m_weight, rep(0,6))
							names(results_m_weight)[(length(results_m_weight)-5):length(results_m_weight)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
							results_m_weight[(length(results_m_weight)-5):length(results_m_weight)] <- unlist(pvalues_ds$results_weight[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
							                                                                                                                  "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
							results_m_weight["STAAR-O"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][p_seq]))
							results_m_weight["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num,6*apc_num+9)]))
							results_m_weight["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
							results_m_weight["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
							results_m_weight["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
							results_m_weight["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
							results_m_weight["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

							results_weight_m <- cbind(results_weight_m, results_m_weight)
							colnames(results_weight_m)[i] <- c(i-1)
						}

						#results_weight1
						for(i in 1:ncol(pvalues_m$results_weight1)){
							results_m_weight1 <- pvalues_m$results_weight1[-c(1,2),i]
							results_m_weight1 <- unlist(pvalues_m$results_weight1[,i][c(5:length(pvalues_m$results_weight1[,i]), 4,3)])
							names(results_m_weight1) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

							results_m_weight1 <- c(results_m_weight1, rep(0,6))
							names(results_m_weight1)[(length(results_m_weight1)-5):length(results_m_weight1)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
							results_m_weight1[(length(results_m_weight1)-5):length(results_m_weight1)] <- unlist(pvalues_ds$results_weight1[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
							                                                                                                                      "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
							results_m_weight1["STAAR-O"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][p_seq]))
							results_m_weight1["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num,6*apc_num+9)]))
							results_m_weight1["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
							results_m_weight1["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
							results_m_weight1["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
							results_m_weight1["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
							results_m_weight1["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

							results_weight1_m <- cbind(results_weight1_m, results_m_weight1)
							colnames(results_weight1_m)[i] <- c(i-1)
						}

						#results_weight2
						for(i in 1:ncol(pvalues_m$results_weight2)){
							results_m_weight2 <- pvalues_m$results_weight2[-c(1,2),i]
							results_m_weight2 <- unlist(pvalues_m$results_weight2[,i][c(5:length(pvalues_m$results_weight2[,i]), 4,3)])
							names(results_m_weight2) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

							results_m_weight2 <- c(results_m_weight2, rep(0,6))
							names(results_m_weight2)[(length(results_m_weight2)-5):length(results_m_weight2)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
							results_m_weight2[(length(results_m_weight2)-5):length(results_m_weight2)] <- unlist(pvalues_ds$results_weight2[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
							                                                                                                                      "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
							results_m_weight2["STAAR-O"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][p_seq]))
							results_m_weight2["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num,6*apc_num+9)]))
							results_m_weight2["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
							results_m_weight2["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
							results_m_weight2["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
							results_m_weight2["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
							results_m_weight2["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

							results_weight2_m <- cbind(results_weight2_m, results_m_weight2)
							colnames(results_weight2_m)[i] <- c(i-1)
						}

						rownames(pvalues_m$weight_all_1) <-rownames(pvalues_m$weight_all_2) <- unique(obj_nullmodel$pop.groups)

						results <- list(results,
						                weight_all_1 = pvalues_m$weight_all_1,
						                weight_all_2 = pvalues_m$weight_all_2,
						                results_weight = results_weight_m,
						                results_weight1 = results_weight1_m,
						                results_weight2 = results_weight2_m)
					}
				}else
				{
					results_m <- c(results[1,],rep(0,2))
					names(results_m)[(length(results_m)-1):length(results_m)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
					results_m[(length(results_m)-1):length(results_m)] <- results[2,c("Burden(1,25)","Burden(1,1)")]

					## check whether the p-values is NA. If so, set NA equals 1.
					if(is.na(results_m[(length(results_m)-1)]))
					{
						results_m[(length(results_m)-1)] <- 1
					}

					if(is.na(results_m[length(results_m)]))
					{
						results_m[length(results_m)] <- 1
					}

					apc_num <- (length(results_m)-10)/2
					p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),(length(results_m)-6):(length(results_m)-5))

					## calculate STAAR-B
					pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
					if(sum(is.na(pvalues_sub))>0)
					{
						if(sum(is.na(pvalues_sub))==length(pvalues_sub))
						{
							results_m["STAAR-B"] <- 1
						}else
						{
							## not all NAs
							pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
							if(sum(pvalues_sub[pvalues_sub<1])>0)
							{
								## not all ones
								results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])

							}else
							{
								results_m["STAAR-B"] <- 1

							}
						}
					}else
					{
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
						}else
						{
							results_m["STAAR-B"] <- 1
						}
					}

					## calculate STAAR-B(1,25)
					pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
					if(sum(is.na(pvalues_sub))>0)
					{
						if(sum(is.na(pvalues_sub))==length(pvalues_sub))
						{
							results_m["STAAR-B(1,25)"] <- 1
						}else
						{
							## not all NAs
							pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
							if(sum(pvalues_sub[pvalues_sub<1])>0)
							{
								## not all ones
								results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])

							}else
							{
								results_m["STAAR-B(1,25)"] <- 1

							}
						}
					}else
					{
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
						}else
						{
							results_m["STAAR-B(1,25)"] <- 1
						}
					}

					## calculate STAAR-B(1,1)
					pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
					if(sum(is.na(pvalues_sub))>0)
					{
						if(sum(is.na(pvalues_sub))==length(pvalues_sub))
						{
							results_m["STAAR-B(1,1)"] <- 1
						}else
						{
							## not all NAs
							pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
							if(sum(pvalues_sub[pvalues_sub<1])>0)
							{
								## not all ones
								results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])

							}else
							{
								results_m["STAAR-B(1,1)"] <- 1

							}
						}
					}else
					{
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
						}else
						{
							results_m["STAAR-B(1,1)"] <- 1
						}
					}

					results_ds <- c()
					results_ds <- rbind(results_ds,results[2,])

					results <- c()
					results <- rbind(results,results_m)
				}
			}
		}
	}else
	{
		results <- c()
		results_ds <- c()
	}

	#Assign results_weight objects across functional categories
	categories <- c("plof", "plof_ds", "ds", "synonymous")
	if(use_ancestry_informed == TRUE & find_weight == TRUE & !use_SPA){
		for(k in categories){
			if(!is.null(eval(as.name(paste0("results_", k))))){
				pvalues <- eval(as.name(paste0("pvalues_", k)))

			results_weight <- results_weight1 <- results_weight2 <- c()
			for(i in 1:ncol(pvalues$results_weight)){
				results_weight_temp <- pvalues$results_weight[-c(1,2),i]
				results_weight_temp <- unlist(pvalues$results_weight[,i][c(5:length(pvalues$results_weight[,i]), 4,3)])
				names(results_weight_temp) <- colnames(eval(as.name(paste0("results_", k))))[-c(1:5)]

				results_weight <- cbind(results_weight, results_weight_temp)
				colnames(results_weight)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight1)){
				results_weight_temp1 <- pvalues$results_weight1[-c(1,2),i]
				results_weight_temp1 <- unlist(pvalues$results_weight1[,i][c(5:length(pvalues$results_weight1[,i]), 4,3)])
				names(results_weight_temp1) <- colnames(eval(as.name(paste0("results_", k))))[-c(1:5)]

				results_weight1 <- cbind(results_weight1, results_weight_temp1)
				colnames(results_weight1)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight2)){
				results_weight_temp2 <- pvalues$results_weight2[-c(1,2),i]
				results_weight_temp2 <- unlist(pvalues$results_weight2[,i][c(5:length(pvalues$results_weight2[,i]), 4,3)])
				names(results_weight_temp2) <- colnames(eval(as.name(paste0("results_", k))))[-c(1:5)]

				results_weight2 <- cbind(results_weight2, results_weight_temp2)
				colnames(results_weight2)[i] <- c(i-1)
			}
			rownames(pvalues$weight_all_1) <- rownames(pvalues$weight_all_2) <- unique(obj_nullmodel$pop.groups)
			assign(paste0("results_", k), list(eval(as.name(paste0("results_", k))),
			                                   weight_all_1 = pvalues$weight_all_1,
			                                   weight_all_2 = pvalues$weight_all_2,
			                                   results_weight = results_weight,
			                                   results_weight1 = results_weight1,
			                                   results_weight2 = results_weight2))
			}
		}
	}

	results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)

	seqResetFilter(genofile)

	return(results_coding)
}

