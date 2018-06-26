## preparation

agi_genes = c("SLC2A1","SLC16A1","H1F1A","LAMC2","COL17A1","ITGB1","AREG","NRG1","EGFR","CDH3","KRT16","KRT17","MCM2","MCM10","CDC7","CDKN2A","E2F2","RPA2","AKR1C1","AKR1C3","ALDH3A1","ICOS","CD8A","LAG3","HLA-DRA","VIM","MMP9")

# save

eset_2 = eset

colnames(eset) = str_replace(colnames(eset),".CEL","")

sig_tab = read.table("~/Dropbox/PhD/Koop_Klinghammer/Misc/Keck_et_al_821_genes.csv",sep="\t",header=T)
sig_tab = apply( sig_tab, MARGIN = 2, FUN = str_replace, ",", ".")
tmp = sig_tab[,1]
sig_tab = apply( sig_tab[,2:4], MARGIN = 2, FUN = as.numeric)
rownames(sig_tab) = tmp

rownames(exprs(eset) ) = mget( rownames(exprs(eset)), hgu133plus2SYMBOL )

### mapping

map_821 = which(rownames(eset) %in% rownames(sig_tab))
eset2 = eset[map_821,]


#map_27 = which(rownames(eset) %in% agi_genes)
#eset3 = eset[map_27,]

### calculation of key values such as mean expression ###

eset2_avg = apply( exprs(eset2), MARGIN=1, FUN = mean  )
#eset3_avg = apply( exprs(eset3), MARGIN=1, FUN = mean  )

myfun = function(vec, agg_list){ vec=aggregate(vec ~ agg_list, FUN = mean)[,2]; return(vec) } 

d2 = apply( exprs(eset2), MARGIN=2, FUN = myfun, rownames( eset2 ) )
rownames(d2) = unique(rownames(eset2))
sig_tab = sig_tab[ which( rownames(sig_tab) %in% rownames(d2) )  ,]

# for 27 genes from keck

d27 = apply( exprs( eset3 ), MARGIN = 2, FUN = myfun, rownames( eset3 )  )
d27_avg = apply( d27, MARGIN=1, FUN = mean  )
d27_minus_avg = d27 - d27_avg

# for 821 genes

means = apply( exprs(eset2), MARGIN = 1, FUN = mean )
d3 = matrix( exprs(eset2) - as.double(means), ncol = dim(exprs(eset2))[2] ) 
d_norm = apply( d3 , MARGIN = 1, FUN = function(vec){ 

  vec = as.double(unlist(vec))
  if( min(vec) < 0 ){  
    vec = vec + abs(min(vec)) 
  } else { 
    vec = vec - abs(min(vec))
  }
  
  vec = vec / max(vec); 
  if (! is.na(max(vec))){
    
    vec = vec - .5; 
    vec = vec * 2; 
    return( vec )  
  }
})
d_norm = matrix( unlist(d_norm), ncol = dim( exprs(eset2) )[2] )


dif_fun = function( vec, compare_mat ){ dif = apply( sqrt( (compare_mat - vec) **2 )  , FUN = sum, MARGIN  = 2 ); return(dif) }
diffs = apply(d3, FUN=dif_fun, MARGIN = 2, sig_tab)

library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
m=colorRampPalette(colors = c("green","black","red"))( 10 )
heatmap.3(d3, col = m)

### Clustering

library("amap")
clusters = Kmeans(x = t(exprs(eset)), centers = 3, method="euclidean", iter.max=100)
table(clusters$cluster)
