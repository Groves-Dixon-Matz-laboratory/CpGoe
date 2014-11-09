Sample_Gene_Expression_Data_README_file

Provides information on the data table 'Sample_Gene_Expression_Data.txt'

The table is organized so that rows hold data on each isogroup (gene)

#Columns
Isogroup: 
	The isogroup (gene). The numbers are associated with the transcriptome
	used to map the RNAseq data. The transcriptome was published by Moya et al. 2012
	and is freely available on the web.
	citation:
	Moya a, Huisman L, Ball EE, Hayward DC, Grasso LC, Chua CM, Woo HN, Gattuso J-P, 
	Forêt S, Miller DJ: Whole transcriptome analysis of the coral Acropora millepora 
	reveals complex responses to CO₂-driven acidification during the initiation of 
	calcification. Mol Ecol 2012, 21:2440–54.
KK1, KK10, KK11 etc...
	These column heading represent the samples. Data in these columns are fold change
	differences in expression generated using DESeq for the indicated sample. Sample
	notation is as follows: The first letter represents the site the sample originated 
	from, either Keppel or Orpheus. The second letter represents the site the sample
	was placed at during the duration of the experiment, again either Keppel or Orpheus.
	The number following the two letters indicates the biological replicate number.
	See Sample_Information_Table.txt for more information on individual samples.
mean_expression:
	This column hold the mean expression in terms of fold change across all samples.
adjp.t:
	This column holds the FDR adjusted p values for the effect of transplantation on 
	gene expression. All p values were generated using DESeq.
adjp.o:
	This column holds the FDR adjusted p values for the effect or site of origin on 
	gene expression.
pval.t:
	This column holds the unadjusted p values for the effect of transplant on gene
	expression.
pval.o:
	This column holds the unadjusted p values for the effect of site of origin on
	gene expression.