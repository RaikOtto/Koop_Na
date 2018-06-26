library("stringr")

balanced.centroid = read.table( "~/Koop_Klinghammer/Misc//balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
balanced.centroid_importance = sort(rowSums(abs(balanced.centroid)), decreasing = T)
balanced.centroid = balanced.centroid[ match(names(balanced.centroid_importance),rownames(balanced.centroid)),]

#pure_data = read.table("~/Koop_Klinghammer/Data/35S.14_03_2018.normalized.tsv", sep ="\t", header = T, row.names = 1)

### Preparation

rownames(pure_data) = str_to_upper(rownames(pure_data))
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "\\_", "")
rownames(pure_data) = str_replace_all(rownames(pure_data), pattern = "-", "")
rownames(balanced.centroid) = str_to_upper(rownames(balanced.centroid))
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "\\_", "")
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "-", "")
colnames(pure_data) = str_replace_all(colnames(pure_data), pattern = "^X", "")

col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)

pure_data = pure_data[row_var > 0, col_var > 0]
dim(pure_data)

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_data = meta_info[match(colnames(pure_data), as.character( meta_info$Name)),]

###

expr = pure_data
source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
table( rownames(expr) %in% rownames(balanced.centroid) )

centroid_genes = rownames(balanced.centroid)
genes_matching_centroid = rownames(expr)[which(  (rownames(expr) %in% rownames(balanced.centroid) ) ) ]
genes_matching_centroid = c( genes_matching_centroid, rep("",length(centroid_genes) - length(genes_matching_centroid)) )
genes_not_matching_centroid = rownames(expr)[which(!(  (rownames(expr) %in% rownames(balanced.centroid) ) )) ]
genes_not_matching_centroid = c(genes_not_matching_centroid, rep("",length(centroid_genes) - length(genes_not_matching_centroid)) )

genes_present = data.frame(
  "Matching_to_centroid" = genes_matching_centroid,
  "Not_matching_to_centroid" = genes_not_matching_centroid,
  "Centroid_genes" = centroid_genes
)
#write.table(genes_present, "~/Koop_Klinghammer/Results/First_results_11S_1C/Centroid_genes_And_not_centroid_genes.tsv",sep ="\t",quote = F, row.names = F)


### centroid classification

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )
expr2bc = centroid2expr( balanced.centroid[,], expr )
colnames(expr2bc$correlation) = c("Sample","Subtype","Correlation","P_value")
class_data = as.data.frame(expr2bc$correlation)

meta_match = match( class_data$Sample, meta_info$Name, nomatch = 0 )
meta_info$Subtype[meta_match] = as.character( class_data$Subtype )
meta_info$P_value[meta_match] = as.double( as.character( class_data$P_value ) )


#write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names =F)

#s_match = match( meta_data$Name, meta_info$Name, nomatch = 0 )
#meta_info$Subtype_Kmeans[s_match] = meta_data$Subtype_Kmeans
