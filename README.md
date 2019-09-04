# metaITH  

metaITH is an R package for analyzing intratumor heterogeneity at the DNA, RNA, and immune levels. Options include a comprehensive signature analysis looking at the anti-PD1, apoptosis, epithelial-mesenchymal, hypoxia, pemetrexed resistance, and proliferation signatures, a dendrogram analysis at the DNA, RNA, and immune levels, and a multiregion divergence and diversity analysis using data at the DNA, RNA, and immune levels. 


# INSTALLATION 

Must have R version 3.5 or higher 
```r
install.packages(devtools) 
library(devtools) 
install_packages("Elise-Merritt/metaITH") 
library(metaITH) 
```

To run metaITH using example datasets, the metaITH package must be downloaded and unzipped. The example datasets are in the folder input_files. 



# USAGE 
  

## Signature Analysis
Functions calculate specific signature scores for the anti-PD1, apoptosis, epithelial-mesenchymal, hypoxia, pemetrexed resistance, and proliferation signatures, that calculate all the above signature scores, and that calculate the signature scores for a geneset supplied by the user. All functions also create a heatmap of the scores.
Functions that calculate signature scores for a specific signature take as input a matrix of z-scores that can be created by running the function `z_score_calculations()`  

Example-calculation of proliferation scores:  
```r
z_score_calculations("RNA_matrix.txt")
proliferation_signature("z-scores_matrix.txt")
read.table("Proliferation_score.txt", sep="\t", header=T)[,1:4]
```

## Dendrogram Analysis
Functions create distance matrices from DNA, RNA, or immune data and use the distance matrices to create unrooted dendrograms, as well as SNV heatmaps.  

Example-creation of DNA dendrogram:  
```r
dna_dendrograms("DNA_dendro_list.txt")
read.table("DNA_distance_matrix_DNA_S1_frequency_matrix.txt", sep="\t", header=T)
snv_heatmaps("DNA_dendro_list.txt")
```

## Multi-region Divergence and Diversity Analysis
Function provides an estimate of the ratio of intra-patient regional diversity to tumor-benign tissue divergence, using the output of the dendrogram analysis functions. 

Example:
```r
multi_level_divergence_diversity("sample_names.txt")
```

# Citation  
A. Sharma, E. Merritt, A. Cruz, C. Jiang, H. Sarkodie, Z. Zhou, J. Malhotra, S. De, Non-genetic intra-tumor heterogeneity is a major predictor of phenotypic heterogeneity and ongoing evolutionary dynamics in lung tumors, Bioxriv (2019).







