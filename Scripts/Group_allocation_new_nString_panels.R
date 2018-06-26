library("pamr")
library("stringr")

balanced.centroid = read.table( "~//Koop_Klinghammer/Misc//balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)

# pure data

#groups = c("CL","CL","BA","MS","BA","BA","BA","MS","MS","BA","MS","BA","MS","BA","BA","BA","CL","BA","BA","BA","CL","CL","CL","BA","BA","BA","CL","CL")

pure_data = read.table( "~//Koop_Klinghammer/upload/Output/Expression.tab", header = T, sep = "\t", stringsAsFactors = F)

hgnc_genes = pure_data$HGNC

colnames( pure_data ) = str_replace( colnames( pure_data ), pattern = "(^X)", ""  )

names_input = str_replace( colnames( pure_data ), pattern = "(^X)", "" )
names_input = str_replace( names_input, pattern = ".CEL", "" )

n_gram_genes = read.table("~/Koop_Klinghammer/Misc/genes_of_interest.tab",sep="\t",header = T)
n_gram_genes = as.character(n_gram_genes$Progression_panel)

n_gram_mat = pure_data[ 
  match( n_gram_genes, hgnc_genes ) ,
]

n_gram_mat = n_gram_mat[,! (colnames(n_gram_mat) %in% c("Probe_id","HGNC","10922","10980A"))]

index_na = apply(n_gram_mat,function(vec){return(sum(is.na(vec)))},MARGIN = 1)
n_gram_mat = n_gram_mat[ index_na == 0,]
n_gram_mat = scale(matrix(as.double(unlist(n_gram_mat)),ncol = 28))

genes = hgnc_genes[match( n_gram_genes, hgnc_genes )]
genes = genes[index_na==0]
rownames(n_gram_mat) = genes

mydata = list(
  x = n_gram_mat,
  y = factor( groups ),
  genenames = rownames(n_gram_mat),
  geneid = as.character(1:nrow(n_gram_mat))
)

fit = pamr.train(mydata)

pamr.predict( fit = fit, newx = n_gram_mat, threshold= 1.0)

centroids = as.data.frame( fit$centroids,ncol=3)
centroids$means = rowMeans( centroids )
centroids = centroids[order( centroids$means, decreasing = T),]

write.table(centroids,"~/Koop_Klinghammer/Misc/Explorative_centroids_HSNC_ncounter_panel.tab",sep="\t",quote =F,row.names=T)

# optics

pdf("~/Koop_Klinghammer/Misc/Explorative_centroids_259_gene_clean_listing.pdf")
  pamr::pamr.plotcen(fit = fit, data = mydata,threshold = c(1.0))
dev.off()

final_centroid_list = pamr::pamr.listgenes(fit = fit, data = mydata, threshold = 1.0, genenames = T )
write.table(
  final_centroid_list,
  "~/Koop_Klinghammer/Misc/Explorative_centroids_259_gene_clean_listing.tab",
  sep="\t",
  quote = F,
  row.names = T
)
