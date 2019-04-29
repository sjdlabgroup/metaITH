# metaITH
A Tool for Analyzing Intratumor Heterogeneity at Multiple Levels
README file written by Anchal Sharma
Following packages should be installed:

install.packages(ggplot2)
install.packages(reshape2)
install.packages(ape)
install.packages(phylobase)
install.packages(dplyr)

All the input files should be in the same directory as metaITH.R script.





============
# PHYLOGRAMS
============
Input files:
	DNA_phylo_list.txt - name of sample files. Each file contains a matrix with variant allele frequency of all the variations in each tumor region and normal. Example file (S1_SNV_frequency_matrix.txt) is provided in /example/phylograms/inputs directory. S1 stands for the name of sample. In this example S1 has 3 tumor regions and 1 normal region profiled. Header of the first column should be coords.
	RNA_phylo_list.txt - name of sample files. Each file contains a matrix with expression values (preferably log2(TPM+1)) of all the genes in each tumor region and normal. Example file (S1_RNA_expression_matrix.txt) is provided in /example/phylograms/inputs directory. S1 stands for the name of sample. Input file name could be any.
	Immune_phylo_list.txt - name of sample files. Each file contains a matrix with proportion of immune cells (inferred by CIBERSORT) of all immune cell types in each tumor region and normal. Example file (S1_Immune_CIBERSORT_matrix.txt) is provided in /example/phylograms/inputs directory. S1 stands for the name of sample. Input file name could be any.


Output files:
	Distance matrices: 3 distance matrices for each sample one for each DNA, RNA and Immune. It contains distance between all the regions and normal. Example file (DNA_distance_matrix_S1_SNV_frequency_matrix.txt) is provided in /example/phylograms/outputs directory. S1 stands for the name of sample.
	Tree topology: 3 tree topology files for each sample one for each DNA, RNA and Immune. Example file (DNA_tree_S1_SNV_frequency_matrix.txt) is provided in /example/phylograms/outputs directory. S1 stands for the name of sample.
	Phylogram: 3 phylograms for each sample one for each DNA, RNA and Immune. Example file (DNA_unrooted_phylogram_S1_SNV_frequency_matrix.tiff) is provided in /example/phylograms/outputs directory. S1 stands for the name of sample.	
	SNV vaf heatmaps: heatmap of variant allele frequency (vaf) of all the variations for a sample in all regions will be made. Example file (SNV_heatmap_DNA_S1_frequency_matrix.txt.tiff) is provided in /example/phylograms/outputs directory. S1 stands for the name of sample.	





============
# SIGNATURES
============
Input files:
	RNA_expression_matrix - one matrix containing expression of all genes across all tumor regions and normal of all samples. File should be named: RNA_expression_matrix.txt. Example file (RNA_expression_matrix.txt) is provided in /example/signatures/inputs directory. 
					* Header of first column should be Gene. 
					* Normal regions should have header ending in ".nr", e.g. S1.nr, S2.nr. 
					
	Geneset files - multiple files with list of genes associated with signatures. For example list of proliferation associated genes in one file with a header Gene. Following genesets are already provided in genesets/ directory:
					* mesenchymal_geneset.txt
					* epithelial_geneset.txt
					* proliferation_geneset.txt
					* apoptosis_geneset.txt
					* anti-PD1_favor_geneset.txt
					* pemetrexed_resistance_geneset.txt
					* hyppoxia_geneset.txt
		
Output files:
	z-score matrix: one matrix containing z-scores for all the genes in all the tumor samples. Example file (z-scores_matrix.txt) is provided in /example/signatures/outputs directory. This matrix can also be used for calculating scores for any other geneset than provided here.
	score files: one score file for each of the above mentioned genesets. This has scores for all the samples. Example file (Proliferation_scores.txt) is provided in /example/signatures/outputs directory.
	heatmap: one heatmap .tiff file for each of the above mentioned genesets. Example file (Proliferation_scores_heatmap.tiff) is provided in /example/signatures/outputs directory.	
	




======================================
# MULTI-LEVEL DIVERGENCE AND DIVERSITY
======================================
Input files:
	sample_names.txt: .txt file with names of all samples separated by new line. Example sample_names.txt is given in example/divs-divg/inputs directory. 
	Distance matrices: Distance matrices (DNA, RNA and Immune) generated above will act as input for this code. 
	
Output files:
	png file: one png file with 8 plots.
					* 3 diversity/divergence plots: DNA v/s RNA, DNA v/s Immune, RNA v/s Immune. Each dot represents one sample.
					* 3 diversity plots: DNA v/s RNA, DNA v/s Immune, RNA v/s Immune. Each dot of same color different regions of same sample.
