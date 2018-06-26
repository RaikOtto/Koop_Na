library("stringr")
library("ggplot2")

s_match = match( as.character( colnames(pure_data)), as.character( meta_info$SampleID), nomatch = 0)
meta_data = meta_info[s_match,]
rownames(meta_data) = meta_data$SampleID

#threshold = 0.05
#meta_data$Sig[ which( meta_data$P_value < threshold ) ] = "Sig"
#meta_data$Sig[ which( meta_data$P_value >= threshold ) ] = "Not_Sig"
#meta_data$Subtype[ which( meta_data$P_value >= threshold ) ] = "Not_sig"

#meta_data$Sig[ (as.character(meta_data$Erhaltungstherapie) != "2" ) ] = "Not_sig"
#meta_data = meta_data[meta_data$Sig == "Sig",]

pure_data = pure_data[,colnames(pure_data) %in% meta_data$SampleID]

#meta_data$Subtype_2[! meta_data$Sig_2] = "Not_sig"
#meta_data$P_value_2 = 1-meta_data$P_value_2

#hc.cut <- factoextra::hcut(
#  scale(t(pure_data)),
#  k = 3,
#  hc_method = "ward.D"
#)
#factoextra::fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
#factoextra::fviz_cluster(hc.cut, ellipse.type = "convex")
#meta_data$Subtype = hc.cut$cluster

###

aka3 = list(
  Organ   = c(
      abdomen = "magenta",
      brain = "cyan",
      Femur = "brown",
      nasal = "lightgreen",
      neck = "brown",
      nodal = "darkgreen",
      pharynx = "purple",
      retrosternal = "blue",
      skin = "gray",
      spleen = "darkred",
      stomach = "pink",
      tbd = "white",
      Testis = "red",
      thyroid = "orange",
      tonsil = "yellow"
  ),
  Subtype = c(
      ABC = "brown",
      GCB = "black",
      healthy = "green",
      nodal = "white",
      PTL = "blue",
      tbd = "pink",
      Test_relapse = "gray"
  ),
  relapse = c(
      healthy = "green",
      no_relapse = "white",
      relapse = "red",
      will_relapse = "orange"
  ),
  relapse_dichotomic = c(
    healthy = "green",
    no = "white",
    yes = "red"
  )
)

### vis

cor_mat = cor(pure_data)
pca = prcomp(t(cor_mat));

b_match = match(rownames(pure_data), rownames(balanced.centroid), nomatch = 0) 
cent_clust = t( balanced.centroid[ b_match,] )
clust_data = pure_data[b_match != 0,]

fit = kmeans(t(scale(clust_data)), centers = cent_clust )
meta_data$Subtype_Kmeans = fit$cluster
meta_data$Subtype_Kmeans[meta_data$Subtype_Kmeans == "1"] = "MS"
meta_data$Subtype_Kmeans[meta_data$Subtype_Kmeans == "2"] = "BA"
meta_data$Subtype_Kmeans[meta_data$Subtype_Kmeans == "3"] = "CL"

ggbiplot::ggbiplot(
  pca,
  obs.scale = 1,
  var.scale = 1, 
  #labels = colnames(pure_data),
  alpha = 1,
  groups = as.character( meta_data$Testis ),
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)
#+ theme(legend.position = "top") + geom_point(aes(colour = meta_data$Subtype_Kmeans, size = meta_data$OS^1.5))

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Organ","Subtype","relapse","relapse_dichotomic","Testis")],
  show_rownames = F,
  show_colnames = T,
  clustering_method = "ward.D2",
  annotation_colors = aka3
)
