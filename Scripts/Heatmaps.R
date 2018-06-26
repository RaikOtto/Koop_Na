
### run group_allocation.R before

### heatmap

library(devtools)
install_github("vqv/ggbiplot")
require(RColorBrewer)

cols = colorRampPalette(
  brewer.pal(
    10,
    "RdBu"
  )
)(256)

subtype_top_bar = as.matrix(
  c(  
    colorRampPalette(
      colors = c("blue")
    )( length( groups ) )
  ) 
)

subtype_top_bar[ groups == "BA",  ] = "yellow"
subtype_top_bar[ groups == "CL",  ] = "green"
#subtype_top_bar[ !( classifier_groups %in% c("CL","BA","MS") ),  ] = "black"

colnames(subtype_top_bar) = c("Subtype")

#cols = c(
#  rep("red", 33),
#  rep("black", 34),
#  rep("green", 33)
#)

pdf("~/Koop_Klinghammer/gene expression paper/Grafiken/Revision/Scaled_row_mean_minus_sample_expression.pdf")
heatmap.3( 
    mean_mat,
    Rowv = T,
    Colv = F,
    dendrogram = "none",
    trace = "none",
    scale = "none",
    col = rev(cols),
    ColSideColors = subtype_top_bar,
    labRow = ""
)

legend("left",      # location of the legend on the heatmap plot
       legend = c("BA", "CL", "MS"), # category labels
       col = c("yellow", "green", "blue"),  # color key
       lty= 1,             # line style
       lwd = 10
)
dev.off()
  
pdf("~/Koop_Klinghammer/gene expression paper/Grafiken/Revision/Absolute_expression_716(821)_genes.pdf")
heatmap.3(
    result,
    Rowv = T,
    Colv = T,
    dendrogram = "none",
    trace = "none",
    labRow = "",
    col = rev(cols),
    ColSideColors = subtype_top_bar
)

legend("left",      # location of the legend on the heatmap plot
       legend = c("BA", "CL", "MS"), # category labels
       col = c("yellow", "green", "blue"),  # color key
       lty= 1,             # line style
       lwd = 10
)
dev.off()

plot_ly(
  z = mean_mat,
  type = "heatmap",
  colors = c("green","lightgreen","black","orange","red")
)

vals <- unique(scales::rescale(c(mean_mat)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = mean_mat, colorscale = colz, type = "heatmap")

## ggplot2

library(ggplot2)

p = ggplot(
  result_mat,
  aes(variable, Name)
)
p = p + geom_tile(aes(fill = rescale) + colour = "white")
p = p + scale_fill_gradient(low = "white", high = "steelblue")