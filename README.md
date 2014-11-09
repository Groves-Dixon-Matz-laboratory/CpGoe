CpGoe
=====

Scripts and instructions for running CpGoe analyses for inferring DNA methylation


Scripts:
getCpGoe.py:
	Extracts CpG count data from a fasta file
DESeq_RT_final.R
	Implements DESeq to measure differential expression due to transplantion and origin
CpG_distribution.R
	Analyzes CpG and expression data and generates figures

Data Files:
Origin_Differential_Expression_Data.txt
	Table of Expression Difference (log2.difference) between colony halves grouped by site of origin
Transplant_Differential_Expresion_Data.txt
	Table of Expression Difference (log2.difference) between colony halves grouped by transplant site
Sample_Gene_Expression_Data.txt
	Table of normalized gene expression for each gene (isogroup) for each sample including global mean expression
Sample_Information_Table.txt
	Table of summary information on each sample