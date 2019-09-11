
############################################################################### Signatures #####################################################################
# Signature Analysis
# =================
#' Runs a signature analysis on an RNA expression matrix
#' 
#' Takes in RNA expression matrix containing the expression level of all genes in all tumor samples and normals.
#' Creates a z-score matrix and uses this matrix to find the hypoxia, proliferation, apoptosis, drug-resistance, emt, and anti-PD1 favor scores for each tumor sample or region. Outputs a table and heatmap of these scores.
#' @usage signature_analysis(rna_expression_matrix_file, type)
#' @param rna_expression_matrix_file File containing expression level of all genes in all tumor samples. Header of first column should be "Gene" and each normal file should end in ".nr"
#' @param type String specifying type of analysis. Options are "all", "hypoxia", "proliferation", "apoptosis", "drug-resistance", "emt", "anti-PD1". Default parameter is "all".
#' @example signature_analysis(RNA_matrix.txt, "all")
#' @export 
signature_analysis=function(rna_expression_matrix_file, type="all"){
orimatrix=read.table(rna_expression_matrix_file, sep="\t", header=T)
nr_matrix=dplyr::select(orimatrix, ends_with("nr"))
orimatrix=dplyr::select(orimatrix, -ends_with("nr"))
names(orimatrix)[1]="Gene"
nr_matrix$Gene=orimatrix$Gene
nr_matrix=nr_matrix[c("Gene", setdiff(names(nr_matrix), "Gene"))]

# For normal file, find mean and standard deviation
Mean=rowMeans(nr_matrix[,2:ncol(nr_matrix)])
SD=rowSds(data.matrix(nr_matrix), cols=2:ncol(nr_matrix))
nr_matrix=cbind(nr_matrix,Mean,SD)
nmatrix=data.frame(nr_matrix$Gene, nr_matrix$Mean, nr_matrix$SD)
names(nmatrix)=c("Gene","Mean", "SD")

# Join normal file calculations to tumor matrix
orimatrix=dplyr::left_join(orimatrix, nmatrix, by="Gene")

# Find z-scores and make new file for output
zscore_matrix=data.frame(matrix(ncol=ncol(orimatrix)-2, nrow=nrow(orimatrix)))
zscore_matrix[,1]=orimatrix$Gene
names(zscore_matrix)[1]="Gene"
for(i in 2:(ncol(orimatrix)-2)){
  sample_name=colnames(orimatrix)[i]
  names(zscore_matrix)[i]=paste0(sample_name, ".")
  zscore_matrix[,i]=(orimatrix[,i]-orimatrix$Mean)/orimatrix$SD
}
write.table(zscore_matrix, "z-scores_matrix.txt", sep="\t", quote=F, row.names=F, col.names=T)

if(type=="all"){
  hypoxia_signature("z-scores_matrix.txt")
  proliferation_signature("z-scores_matrix.txt")
  apoptosis_signature("z-scores_matrix.txt")
  drug_resistance_signature("z-scores_matrix.txt")
  emt_signature("z-scores_matrix.txt")
  anti_pd1_favor_signature("z-scores_matrix.txt")
}else if(type=="hypoxia"){
  hypoxia_signature("z-scores_matrix.txt")
}else if(type=="proliferation"){
  proliferation_signature("z-scores_matrix.txt")
}else if(type=="apoptosis"){
  apoptosis_signature("z-scores_matrix.txt")
}else if(type=="drug-resistance"){
  drug_resistance_signature("z-scores_matrix.txt")
}else if(type=="emt_signature"){
  emt_signature("z-scores_matrix.txt")
}else if(type=="anti-PD1"){
  anti_pd1_favor_signature("z-scores_matrix.txt")
}else{
  hypoxia_signature("z-scores_matrix.txt")
  proliferation_signature("z-scores_matrix.txt")
  apoptosis_signature("z-scores_matrix.txt")
  drug_resistance_signature("z-scores_matrix.txt")
  emt_signature("z-scores_matrix.txt")
  anti_pd1_favor_signature("z-scores_matrix.txt")
}
}



