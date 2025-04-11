#' Individual-variant analysis using score test
#'
#' The \code{Individual_Analysis} function takes in chromosome, starting location, ending location, an user-defined variant list for
#' ancestry-informed analyses, the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype (including imbalanced case-control design) and each individual variant in a genetic region by using score test.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait score test p-values by leveraging
#' the correlation structure between multiple phenotypes.
#' For ancestry-informed analysis, the results correspond to ensemble p-values across base tests,
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param start_loc starting location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param end_loc ending location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param individual_results the data frame of (significant) individual variants of interest for ancestry-informed analysis.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param mac_cutoff the cutoff of minimum minor allele count in
#' defining individual variants (default = 20).
#' @param subset_variants_num the number of variants to run per subset for each time (default = 5e3).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint approximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integer specifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05)
#' @param use_ancestry_informed logical: is ancestry-informed association analysis used to estimate p-values (default = FALSE).
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @return A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region, or as provided in \code{individual_results}
#' for ancestry-informed variant analysis. The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' If \code{find_weight} is TRUE, returns a list containing the ancestry-informed score test p-values and the estimated effect size of the minor allele for each individual variant provided in \code{individual_results}.
#' The ensemble weights under two sampling scenarios and p-values under scenarios 1, 2, and combined for each base test are saved as well.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

Individual_Analysis <- function(chr,start_loc=NULL,end_loc=NULL,individual_results=NULL,genofile,obj_nullmodel,mac_cutoff=20,subset_variants_num=5e3,
                                QC_label="annotation/filter",variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                tol=.Machine$double.eps^0.25,max_iter=1000,SPA_p_filter=TRUE,p_filter_cutoff=0.05,
                                use_ancestry_informed=FALSE,find_weight=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)
	n_pheno <- obj_nullmodel$n.pheno

	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}

	## residuals and cov
	residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
	if(SPA_p_filter)
	{
		### dense GRM
		if(!obj_nullmodel$sparse_kins)
		{
			P <- obj_nullmodel$P
		}

		### sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Sigma_i <- obj_nullmodel$Sigma_i
			Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
			cov <- obj_nullmodel$cov
		}
	}

	## SPA
	if(use_SPA)
	{
		muhat <- obj_nullmodel$fitted.values

		if(obj_nullmodel$relatedness)
		{
			if(!obj_nullmodel$sparse_kins)
			{
				XW <- obj_nullmodel$XW
				XXWX_inv <- obj_nullmodel$XXWX_inv
			}else
			{
				XW <- as.matrix(obj_nullmodel$XSigma_i)
				XXWX_inv <- as.matrix(obj_nullmodel$XXSigma_iX_inv)
			}
		}else
		{
			XW <- obj_nullmodel$XW
			XXWX_inv <- obj_nullmodel$XXWX_inv
		}
	}else
	{
		### dense GRM
		if(!obj_nullmodel$sparse_kins)
		{
			P <- obj_nullmodel$P
		}

		### sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Sigma_i <- obj_nullmodel$Sigma_i
			Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
			cov <- obj_nullmodel$cov
		}
	}

	## get SNV id
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

	results <- c()

	if(use_ancestry_informed)
	{
		results <- AI_Individual_Analysis(chr=chr,individual_results=individual_results,genofile=genofile,
		                                  obj_nullmodel=obj_nullmodel,QC_label=QC_label,variant_type=variant_type,
		                                  geno_missing_imputation=geno_missing_imputation,find_weight=find_weight)
		return(results)
	}

	position <- as.numeric(seqGetData(genofile, "position"))

	variant.id <- seqGetData(genofile, "variant.id")
	is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
	SNV.id <- variant.id[is.in]
	
	results <- c()
	
	if(length(SNV.id) == 0)
	{
	  return(results)
	}
	
	## get AF, AC
	seqSetFilter(genofile, variant.id = SNV.id, sample.id = phenotype.id)
	AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
	REF_AF <- AF_AC_Missing$af
	REF_AC <- AF_AC_Missing$ac
	Missing_rate <- AF_AC_Missing$miss
	ALT_AC <- 2*round(samplesize*(1-Missing_rate))-REF_AC
	MAC <- ifelse(REF_AC>=ALT_AC,ALT_AC,REF_AC)
	
	is.include <- !((MAC<mac_cutoff) | is.na(MAC))
	SNV.id <- SNV.id[is.include]
	REF_AF <- REF_AF[is.include]
	Missing_rate <- Missing_rate[is.include]
	rm(AF_AC_Missing,is.include)
	gc()
	
	seqResetFilter(genofile)

	subset.num <- ceiling(length(SNV.id)/subset_variants_num)

	results <- c()
	
	if(subset.num == 0)
	{
		return(results)
	}

	for(kk in 1:subset.num)
	{
		if(kk < subset.num)
		{
			is.in <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
		}
		if(kk == subset.num)
		{
			is.in <- ((kk-1)*subset_variants_num+1):length(SNV.id)
		}

	  ## Genotype imputation and flip
	  REF_AF.in <- REF_AF[is.in]
	  Missing_rate.in <- Missing_rate[is.in]
	  Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=SNV.id[is.in],
	                                        sample.id=phenotype.id,
	                                        REF_AF=REF_AF.in,Missing_rate=Missing_rate.in)
	  Geno <- Genotype_sp$Geno
	  results_information <- Genotype_sp$results_information
	  MAF.in <- results_information$MAF
	  ALT_AF.in <- results_information$ALT_AF
	  Missing_rate.in <- results_information$Missing_rate
	  MAC.in <- round(2*MAF.in*(1-Missing_rate.in)*samplesize)
	  
	  CHR <- results_information$CHR
	  position <- results_information$position
	  REF <- results_information$REF
	  ALT <- results_information$ALT
	  N <- rep(samplesize,length(CHR))
	  
	  if (geno_missing_imputation == "mean")
	  {
	    Geno <- na.replace.sp(Geno,m=2*MAF.in)
	  }
	  if (geno_missing_imputation == "minor")
	  {
	    Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
	    MAF.in <- MAC.in/(2*samplesize)
	  }

		if(!all(CHR==chr))
		{
			warning("chr does not match the chromosome of genofile (the opened aGDS)!")
		}

		if((use_SPA)&!SPA_p_filter)
		{
			if(length(MAF.in)>=1)
			{
			  Geno <- as.matrix(Geno)

				pvalue <- Individual_Score_Test_SPA(Geno,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)

				results_temp <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,ALT_AF=ALT_AF.in,
				                           MAF=MAF.in,N=N,pvalue=pvalue)

				results <- rbind(results,results_temp)
			}
		}else
		{
		  ## Common_variants or variants with relatively high missing rate
		  is.common_highmissing <- (MAF.in>=0.01) | (Missing_rate.in>=0.01)
		  if(sum(is.common_highmissing)>=1)
		  {
		    Geno_common <- Geno[,is.common_highmissing,drop=FALSE]
		    
		    CHR_common <- CHR[is.common_highmissing]
		    position_common <- position[is.common_highmissing]
		    REF_common <- REF[is.common_highmissing]
		    ALT_common <- ALT[is.common_highmissing]
		    MAF_common <- MAF.in[is.common_highmissing]
		    ALT_AF_common <- ALT_AF.in[is.common_highmissing]
		    N_common <- N[is.common_highmissing]
		    
		    ## Split into small chunks to run
		    subset_variants_num_common <- 200
		    subset.num_common <- ceiling(length(CHR_common)/subset_variants_num_common)
		    
		    for(kk_common in 1:subset.num_common)
		    {
		      if(kk_common < subset.num_common)
		      {
		        is.in_common_subset <- ((kk_common-1)*subset_variants_num_common+1):(kk_common*subset_variants_num_common)
		      }
		      if(kk_common == subset.num_common)
		      {
		        is.in_common_subset <- ((kk_common-1)*subset_variants_num_common+1):length(CHR_common)
		      }
		      
		      Geno_common_subset <- Geno_common[,is.in_common_subset,drop=FALSE]
		      
		      CHR_common_subset <- CHR_common[is.in_common_subset]
		      position_common_subset <- position_common[is.in_common_subset]
		      REF_common_subset <- REF_common[is.in_common_subset]
		      ALT_common_subset <- ALT_common[is.in_common_subset]
		      MAF_common_subset <- MAF_common[is.in_common_subset]
		      ALT_AF_common_subset <- ALT_AF_common[is.in_common_subset]
		      N_common_subset <- N_common[is.in_common_subset]
		      
		      ## sparse GRM
		      if(obj_nullmodel$sparse_kins)
		      {
		        if(n_pheno == 1)
		        {
		          Score_test <- Individual_Score_Test_sp(Geno_common_subset, Sigma_i, Sigma_iX, cov, residuals.phenotype)
		        }else
		        {
		          Geno_common_subset <- Diagonal(n = n_pheno) %x% Geno_common_subset
		          Score_test <- Individual_Score_Test_sp_multi(Geno_common_subset, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
		        }
		      }
		      
		      ## dense GRM
		      if(!obj_nullmodel$sparse_kins)
		      {
		        if(n_pheno == 1)
		        {
		          Score_test <- Individual_Score_Test_sp_denseGRM(Geno_common_subset, P, residuals.phenotype)
		        }else
		        {
		          Geno_common_subset <- Diagonal(n = n_pheno) %x% Geno_common_subset
		          Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_common_subset, P, residuals.phenotype, n_pheno)
		        }
		      }
		      
		      ## SPA approximation for small p-values
		      if(use_SPA)
		      {
		        pvalue <- exp(-Score_test$pvalue_log)
		        
		        if(sum(pvalue < p_filter_cutoff)>=1)
		        {
		          is.common_subset_SPA <- as.vector(pvalue < p_filter_cutoff)
		          Geno_common_subset_SPA <- Geno_common_subset[,is.common_subset_SPA,drop=FALSE]
		          Geno_common_subset_SPA <- as.matrix(Geno_common_subset_SPA)
		          
		          pvalue_SPA <- Individual_Score_Test_SPA(Geno_common_subset_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
		          
		          pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
		        }
		      }
		      
		      if(use_SPA)
		      {
		        results_temp <- data.frame(CHR=CHR_common_subset,POS=position_common_subset,REF=REF_common_subset,ALT=ALT_common_subset,ALT_AF=ALT_AF_common_subset,
		                                   MAF=MAF_common_subset,N=N_common_subset,pvalue=pvalue)
		      }else
		      {
		        if(n_pheno == 1)
		        {
		          results_temp <- data.frame(CHR=CHR_common_subset,POS=position_common_subset,REF=REF_common_subset,ALT=ALT_common_subset,ALT_AF=ALT_AF_common_subset,MAF=MAF_common_subset,N=N_common_subset,
		                                     pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
		                                     Score=Score_test$Score,Score_se=Score_test$Score_se,
		                                     Est=Score_test$Est,Est_se=Score_test$Est_se)
		        }else
		        {
		          results_temp <- data.frame(CHR=CHR_common_subset,POS=position_common_subset,REF=REF_common_subset,ALT=ALT_common_subset,ALT_AF=ALT_AF_common_subset,MAF=MAF_common_subset,N=N_common_subset,
		                                     pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
		          results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
		          colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
		        }
		      }
		      results <- rbind(results,results_temp)
		    }
		  }


		  ## Rare_variants with relatively low missing rate
		  is.rare_lowmissing <- (MAF.in<0.01) & (Missing_rate.in<0.01)
			if(sum(is.rare_lowmissing)>=1)
			{
			  Geno_rare <- Geno[,is.rare_lowmissing,drop=FALSE]
			  
			  CHR_rare <- CHR[is.rare_lowmissing]
			  position_rare <- position[is.rare_lowmissing]
			  REF_rare <- REF[is.rare_lowmissing]
			  ALT_rare <- ALT[is.rare_lowmissing]
			  MAF_rare <- MAF.in[is.rare_lowmissing]
			  ALT_AF_rare <- ALT_AF.in[is.rare_lowmissing]
			  N_rare <- N[is.rare_lowmissing]

				## sparse GRM
				if(obj_nullmodel$sparse_kins)
				{
				  if(n_pheno == 1)
				  {
				    Score_test <- Individual_Score_Test_sp(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
				  }else
				  {
				    Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
				    Score_test <- Individual_Score_Test_sp_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
				  }
				}

				## dense GRM
				if(!obj_nullmodel$sparse_kins)
				{
				  if(n_pheno == 1)
				  {
				    Score_test <- Individual_Score_Test_sp_denseGRM(Geno_rare, P, residuals.phenotype)
				  }
				  else
				  {
				    Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
				    Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno)
				  }
				}

				## SPA approximation for small p-values
				if(use_SPA)
				{
				  pvalue <- exp(-Score_test$pvalue_log)
				  
				  is.rare_SPA <- as.vector(pvalue < p_filter_cutoff)
				  if(sum(is.rare_SPA)>=1)
				  {
				    pvalue_SPA <- c()
				    Geno_rare_SPA <- Geno_rare[,is.rare_SPA,drop=FALSE]
				    
				    ## Split into small chunks to run
				    subset_variants_num_SPA <- 50
				    subset.num_SPA <- ceiling(sum(is.rare_SPA)/subset_variants_num_SPA)
				    
				    for(kk_SPA in 1:subset.num_SPA)
				    {
				      if(kk_SPA < subset.num_SPA)
				      {
				        is.in_SPA_subset <- ((kk_SPA-1)*subset_variants_num_SPA+1):(kk_SPA*subset_variants_num_SPA)
				      }
				      if(kk_SPA == subset.num_SPA)
				      {
				        is.in_SPA_subset <- ((kk_SPA-1)*subset_variants_num_SPA+1):sum(is.rare_SPA)
				      }
				      
				      Geno_rare_SPA_subset <- Geno_rare_SPA[,is.in_SPA_subset,drop=FALSE]
				      Geno_rare_SPA_subset <- as.matrix(Geno_rare_SPA_subset)
				      
				      if(length(is.in_SPA_subset)>=1)
				      {
				        pvalue_SPA_subset <- Individual_Score_Test_SPA(Geno_rare_SPA_subset,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
				        pvalue_SPA <- c(pvalue_SPA,pvalue_SPA_subset)
				      }
				    }
				    pvalue[is.rare_SPA] <- pvalue_SPA
				  }
				}

				if(use_SPA)
				{
					results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
					                           pvalue=pvalue)

				}else
				{
					if(n_pheno == 1)
					{
						results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
						                           Score=Score_test$Score,Score_se=Score_test$Score_se,
						                           Est=Score_test$Est,Est_se=Score_test$Est_se)
					}
					else
					{
						results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
						results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
						colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
					}
				}

				results <- rbind(results,results_temp)
			}
		}
		seqResetFilter(genofile)
	}

	if(!is.null(results))
	{
		results <- results[order(results[,2]),]
	}

	return(results)
}

