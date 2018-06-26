
### pca

eset_pca = prcomp(
  t( pure_data ),
  center = TRUE,
  scale. = TRUE
)

groups = as.character( expr2bc$cluster )

library("ggbiplot")
g = ggbiplot( 
  eset_pca, 
  obs.scale = 1, 
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  #groups = as.character( cohort_t$Subtype[match_names] ),
  groups = classifier_groups,
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F,
  labels = colnames(pure_data)
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
#pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/final/Classification.pdf")
# print(g)
#dev.off()

# heatmap 3
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
require(RColorBrewer)

cols = colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor   = function(x) as.dist( 1 - cor( t( x ) ) )
hclustAvg = function(x) hclust(  x, method="average" )

#pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/final/Classification.pdf")


#dev.off()

# top bar start

subtype_top_bar = as.matrix(
  c(  
    colorRampPalette(
      colors = c("blue")
    )( length( groups ) )
  ) 
)

subtype_top_bar[ groups == "BA",  ] = "yellow"
subtype_top_bar[ groups == "CL",  ] = "green"
subtype_top_bar[ !( groups %in% c("CL","BA","MS") ),  ] = "black"
#subtype_top_bar = subtype_top_bar[ groups %in% c("","CL","NA"),  ] = "yellow"

colnames(subtype_top_bar) = c("Subtype")

# top bar end

avg_t_plot = avg_t[, !is.na(groups) & (groups != "")]
groups_pca = groups[ !is.na(groups) & (groups != "") ]

rownames(avg_t_plot) = 1:dim(avg_t_plot)[1]

avg_t_pca = prcomp(
  t(avg_t_plot),
  center = TRUE,
  scale. = TRUE
)

library(ggbiplot)
g = ggbiplot(
  avg_t_pca, 
  obs.scale = 1, 
  var.scale = 1, 
  groups = as.character( groups_pca ),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = colnames(avg_t_plot)
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
