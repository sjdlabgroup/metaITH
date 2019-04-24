
############################################################################### Signatures #####################################################################
# z-score calculations
# =================
#' Runs a signature analysis on an RNA expression matrix
#' 
#' Takes in RNA expression matrix containing the expression level of all genes in all tumor samples and normals.
#' Creates a z-score matrix and uses this matrix to find the hypoxia, proliferation, apoptosis, drug-resistance, emt, and anti-PD1 favor scores for each tumor sample or region. Outputs a table and heatmap of these scores.
#' @usage signature_analysis(rna_expression_matrix_file)
#' @param rna_expression_matrix_file File containing expression level of all genes in all tumor samples. Header of first column should be "Gene" and each normal file should end in ".nr"
#' @example signature_analysis(S1_RNA_expresion_matrix)
#' @export 
signature_analysis=function(rna_expression_matrix_file){
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




# Hypoxia signature
# =================
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
hyp <- metaITH:::hyppoxia_geneset 
y<-inner_join(x,hyp, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="Hypoxia score", t(z))
write.table(z2, "Hypoxia_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("Hypoxia_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("Hypoxia_scores_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "darkorchid4") + 
  coord_equal()
dev.off()




# Proliferation signature
# =======================
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
prol <- metaITH:::proliferation_geneset 
y<-inner_join(x,prol, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="Proliferation score", t(z))
write.table(z2, "Proliferation_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("Proliferation_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("Proliferation_scores_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "deeppink3") + 
  coord_equal()
dev.off()




# Apoptosis signature
# ===================
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
apop <- metaITH:::apoptosis_geneset
y<-inner_join(x,apop, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="Apoptosis score", t(z))
write.table(z2, "Apoptosis_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("Apoptosis_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("Apoptosis_scores_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "deeppink3") + 
  coord_equal()
dev.off()




# Drug resistance signature
# =========================
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
resist <- metaITH:::premetrexed_resistance_geneset
y<-inner_join(x,resist, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="Pemetrexed resistance score", t(z))
write.table(z2, "Pemetrexed_resistance_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("Pemetrexed_resistance_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("Pemetrexed_resistance_scores_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "cyan3") + 
  coord_equal()
dev.off()




# EMT signature
# =============
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
epi <- metaITH:::epithelial_geneset 
y<-inner_join(x,epi, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="Epithelial score", t(z))

mesen <- metaITH:::mesenchymal_geneset
a<-inner_join(x,mesen, by="Gene")
b<-colMeans(a[,2:ncol(a)])
b2<-cbind(Gene="Mesenchymal score", t(b))

q <- rbind(z2,b2)
ME <- b-z
ME2 <- cbind(Gene="M-E score", t(ME))

write.table(q, "Epithelial_Mesenchymal_scores.txt", quote=F, sep="\t", row.names=F)
write.table(ME2, "Mesenchymal-Epithelial_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("Mesenchymal-Epithelial_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("Mesenchymal-Epithelial_score_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "yellow3") + 
  coord_equal()
dev.off()



# anti-PD1 favor signature
# ========================
x <- read.table("z-scores_matrix.txt", header=T, sep="\t")
pd1 <- metaITH:::anti_PD1_favor_geneset
y<-inner_join(x,pd1, by="Gene")
z<-colMeans(y[,2:ncol(y)])
z2<-cbind(Gene="anti-PD1 favor score", t(z))
write.table(z2, "anti-PD1_favor_score.txt", quote=F, sep="\t", row.names=F)

x<-read.table("anti-PD1_favor_score.txt", header=T, sep="\t")
x1<-melt(x,id="Gene")
tiff("anti-PD1_favor_scores_heatmap.tiff", units="in", height = 8, width = 8, res=300)
ggplot(x1, aes(variable, Gene, fill=value))+ 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 10),axis.text.y = element_text(size =10,color ="black"), plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom") + 
  scale_fill_gradient(low = "white", high = "darkgreen") + 
  coord_equal()
dev.off()

}



