#' Gene-centric conditional analysis of long noncoding RNA (ncRNA) category using STAAR procedure for imbalance case-control setting
#'
#' The \code{ncRNA_cond_spa} function takes in chromosome, gene name,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between an
#' imbalanced case-control phenotype and the exonic and splicing category of an ncRNA gene by using STAAR procedure.
#' For each ncRNA category, the conditional STAAR-B p-value is a p-value from an omnibus test
#' that aggregated conditional Burden(1,25) and Burden(1,1),
#' together with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @param chr chromosome.
#' @param gene_name name of the ncRNA gene to be analyzed using STAAR procedure.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param rv_num_cutoff_max_prefilter the cutoff of maximum number of variants
#' before extracting the genotype matrix (default = 1e+09).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return A data frame containing the STAAR p-values (including STAAR-O) corresponding to the exonic and splicing category of the given ncRNA gene.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

ncRNA_cond_spa <- function(chr,gene_name,genofile,obj_nullmodel,known_loci,
                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                           rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                           QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                           SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	### known SNV Info
	known_loci_chr <- known_loci[known_loci[,1]==chr,,drop=FALSE]
	known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),,drop=FALSE]


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

	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	## ncRNA SNVs
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

	variant.id.ncRNA <- variant.id[is.in]

	rm(GENCODE.Category)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)

	rm(variant.id.ncRNA)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
	Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

	Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
	Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
	Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

	rm(GENCODE.Info)
	gc()

	rm(GENCODE.Info.split)
	gc()

	variant.id.ncRNA <- seqGetData(genofile, "variant.id")

	seqResetFilter(genofile)

	### Gene
	is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
	is.in <- union(is.in,which(Gene_list_3==gene_name))

	variant.is.in <- variant.id.ncRNA[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- as.character(seqGetData(genofile,"sample.id"))
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- NULL
	if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
	{
		Geno <- seqGetData(genofile, "$dosage")
		Geno <- Geno[id.genotype.match,,drop=FALSE]
	}

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}


	## Genotype Info
	REF_region <- as.character(seqGetData(genofile, "$ref"))
	ALT_region <- as.character(seqGetData(genofile, "$alt"))

	position_region <- as.numeric(seqGetData(genofile, "position"))

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

	## Exclude RV in the region which needed to be adjusted
	if(dim(known_loci_chr)[1]>=1)
	{
		id_exclude <- c()
		for(i in 1:dim(known_loci_chr)[1])
		{
			id_exclude <- c(id_exclude,which((position_region==known_loci_chr[i,2])&(REF_region==known_loci_chr[i,3])&(ALT_region==known_loci_chr[i,4])))
		}

		if(length(id_exclude)>0)
		{
			Geno <- Geno[,-id_exclude]
			Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[-id_exclude,]
		}
	}

	pvalues <- 0
	try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)

	results <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- rep(NA,4)
		results_temp[3] <- "ncRNA_cond"
		results_temp[2] <- chr
		results_temp[1] <- as.character(gene_name)
		results_temp[4] <- pvalues$num_variant

		results_temp <- c(results_temp,pvalues$cMAC,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)

		results <- rbind(results,results_temp)
	}

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
		colnames(results)[dim(results)[2]] <- c("STAAR-B")
	}

	seqResetFilter(genofile)

	return(results)
}

