
# HMPAT for GWAS in which the individual genotype and phenotype can be accessed
# Notes: This function is for HMPAT which improves the power of transcriptome wide association studies by 
#         aggregating multiple expression prediction models
#
# Importantly, data check, including quality control, missing imputation, and effect allele alignement,
#         needs to be done before the analysis
# Version 1.0
#:::INPUT:::
# y (N X 1) is the outcome phenotype in the GWAS
# G2 (N X m) is the matrix for the genotypes of a set of cis-SNPs of a given gene for examining
#         its association with the phenotype
# weight (n X m) is the matrix for weights for these cis-SNPs, which are obtained with various expression prediction models
#         from an external transcriptome panel
# outcome indicates the type of the genotypes and can take values of "C" for continuous phenotypes
#         and "B" for binary phenotypes

library(harmonicmeanp)
HMPAT_individual <- function(y,G2,weight,outcome){
	y <- y
	G2 <- G2
	weight <- weight
	m <- dim(G2)[2]
	K <- dim(weight)[2]
	p_TWAS <- rep(NA,K)
	X2 <- G2%*%weight

	if (outcome== "C") {
	for (j in 1:K)
	{
	p_TWAS[j] <- summary(glm(y~X2[,j],family = gaussian))$coef[2,4]
	}
	}

	if (outcome == "B") {
	for (j in 1:K)
	{
	p_TWAS[j] <- summary(glm(y~X2[,j],family = binomial))$coef[2,4]
	}
	}
	p_HMPAT = as.vector(c(p.hmp(p_TWAS,L=length(p_TWAS))))
	return(list(p_HMPAT=p_HMPAT,p_TWAS=p_TWAS))
}

#HMPAT for GWAS in which only summary statistics can be accessed
# :::INPUT:::
# Z (m X 1) is the z score for each cis-SNP in the GWAS
# G (N X m) is the matrix for the genotypes of reference panel to compute the LD matrix for these cis-SNPs
# weight (n X m) is the matrix for weights for these cis-SNPs, which are obtained with various expression prediction models
#         from an external transcriptome panel
# lambda indicates the value used to compute the shrunk LD matrix

HMPAT_summary <- function(Z,G,weight,lambda){
	Z <- as.vector(Z)
	G <- as.matrix(G)
	weight <- as.matrix(weight)
	K <- dim(weight)[2]
	p_TWAS <- rep(NA,K)
	R <- as.matrix(cor(G))
	R <- lambda * R + diag(1-lambda, nrow=dim(R)[2], ncol=dim(R)[2])
	for (j in 1:K)
	{
	w <- as.matrix(weight[,j])
	z <- t(w)%*%Z/sqrt(t(w)%*%R%*%w)
	p_TWAS[j] <- pnorm(-abs(z))*2
	}
	p_HMPAT = as.vector(c(p.hmp(p_TWAS,L=length(p_TWAS))))
	return(list(p_HMPAT=p_HMPAT,p_TWAS=p_TWAS))
}
