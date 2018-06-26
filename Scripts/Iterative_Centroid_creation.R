library("pamr")

count_data = pure_data[,meta_data$Sig == "Sig"]
colnames(count_data) = str_replace(colnames(count_data),pattern = "^X","")
colnames(count_data) = str_replace(colnames(count_data),pattern = "\\.","_")

count_match = match( colnames(count_data), meta_data$Name, nomatch = 0)

count_data[ count_data == 0] = 1
count_data = log(count_data)

# 6264

centroid2expr = function( centroid, vd ){
  
  gene.sig = intersect( rownames( centroid ), rownames( vd ) )
  vd = t( scale( t( vd[ gene.sig, ] ) ) )
  centroid = centroid[ gene.sig, ]
  vclass = c()
  vcor = c()
  
  for( i in 1:ncol( vd ) ){
    
    d = vd[,i]
    c.cor = c()
    pv = c()
    
    for( j in colnames( centroid ) ){
      
      centroidj = centroid[ , j ]
      corj = cor.test( centroidj, d, use = "complete", method = "pearson" )
      c.cor[ j ] = corj$estimate
      pv[j]=corj$p.value
    }
    
    maxk = which.max(c.cor)
    pub_cor <<- rbind(pub_cor,c.cor)
    group = names( maxk )
    vcor = rbind( vcor, c( colnames( vd )[ i ], group, c.cor[ maxk ], pv[ maxk ] ) )
    
    if( pv[ maxk ] < .05 ){
      vclass[colnames(vd)[i]]=group
    }
  }
  
  
  return( list( overlap.gene = gene.sig, cluster = vclass, correlation = vcor ) )
}
ground_truth_vec <<- as.character(meta_data$Subtype[ match(colnames(count_data), meta_data$Name ) ])
names(ground_truth_vec) = colnames(count_data)

### start

nr_samples_old <<- length(ground_truth_vec)
nr_samples_new <<- length(ground_truth_vec) - 1
iteration <<- 0

while ( nr_samples_old != nr_samples_new ){

  print(iteration)
  nr_samples_old <<- length(ground_truth_vec)
  
  c_match = match(names(ground_truth_vec),colnames(count_data) , nomatch = 0)
  count_data_centroid = count_data[,c_match]
  
  row_var = apply(
    as.matrix(count_data_centroid),
    FUN = function(vec){
      return (var(as.double(vec) ))
    },
    MARGIN = 1
  )
  
  col_var = apply(
    as.matrix(count_data_centroid),
    FUN = function(vec){return (var(vec) )},
    MARGIN = 2
  )
  
  count_data_centroid = as.matrix(count_data_centroid[row_var > 0,col_var > 0])
  print(ncol(count_data_centroid))
  
  ground_truth_vec <<- ground_truth_vec[ names(ground_truth_vec) %in% colnames(count_data_centroid)]
  
  my_data = list(
    x = scale(count_data_centroid),
    y = as.character(ground_truth_vec),
    genenames = rownames(count_data)#,
    #geneid = rownames(count_data)
  )

  fit = pamr.train(
    my_data,
    n.threshold = 30
  )
  min_threshold = tail(which( fit$errors == min(fit$errors)),1)

  centroids <<- as.vector(rep("",ncol(count_data_centroid)))
  centroids = as.data.frame( fit$centroids,ncol=3)
  centroids$means = rowMeans( centroids )
  centroids = centroids[order( centroids$means, decreasing = T),]
  centroids = centroids[,-ncol(centroids)]

  # optics

  #pdf("~/MAPTor_Net_RNA_data/Results/8685/Centroids_9651_IT9.pdf")
    #pamr::pamr.plotcen(fit = fit, data = my_data,threshold = min_threshold)
  #dev.off()

  #{ sink("/dev/null")
  final_centroid_list = pamr::pamr.listgenes(
    fit = fit,
    data = my_data,
    threshold = min_threshold,
    genenames = T
  )
  colnames(final_centroid_list) = str_replace_all(colnames(final_centroid_list) , pattern = "-score", "")

  centroid_hgnc_symbols = rownames()
 #sink(); }
  
  final_centroid_list[,1]=centroid_hgnc_symbols
  final_centroid_list = final_centroid_list[,colnames(final_centroid_list) != "avg_importance"]
  
  write.table(final_centroid_list,
    paste0(sep = "", collapse = "", c("~/MAPTor_NET//Results/Classification/Centroids_IT",iteration,".tsv") ),
  sep="\t",quote =F,row.names=F
  )

  ### classification

  final_centroid_list = read.table(
    paste0(sep = "", collapse = "", c("~/MAPTor_NET//Results/Classification/Centroids_IT",iteration,".tsv") ),
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  final_centroid_list = final_centroid_list[,colnames(final_centroid_list) != "avg_importance"]
  balanced.centroid = as.data.frame(final_centroid_list)[,c(-1,-2)]
  rownames(balanced.centroid) = as.character(final_centroid_list[,2])
  en_r = rownames(balanced.centroid)
  type_c = colnames(balanced.centroid)
  balanced.centroid = t(apply( balanced.centroid, FUN = function(vec){return((as.double(unlist(vec))))} , MARGIN = 1 ))
  balanced.centroid = as.data.frame(balanced.centroid)
  colnames(balanced.centroid) = type_c
  rownames(balanced.centroid) = en_r
  
  # filter for unwanted samples
  
  pub_cor <<- matrix( as.double(), ncol = ncol( balanced.centroid ) )
  avg_val = as.double(rowMeans(abs(balanced.centroid)))
  balanced.centroid = balanced.centroid[order(avg_val, decreasing = T),]
  
  avg_importance = rowSums( abs( balanced.centroid) )
  balanced.centroid = balanced.centroid[order(avg_importance, decreasing = T),]
  balanced.centroid = balanced.centroid[1:100,]
  
  expr2bc = centroid2expr( balanced.centroid, count_data_centroid )
  groups = as.character( expr2bc$correlation[,2] )
  groups = str_replace(groups, pattern = "-score", "")
  colnames(expr2bc$correlation) = c("ID","Group","Cor","PValue")
  
  p_value = as.double(as.character(unlist(expr2bc$correlation[,4])))
  q_value = p.adjust(p_value,"BH")
  correlation = expr2bc$correlation[,3]
  expr2bc$correlation = cbind( expr2bc$correlation , q_value)
  expr2bc$correlation[,4] = p_value
  
  cent_match = match(colnames(count_data_centroid), meta_data$Name, nomatch = 0)
  
  classification_t    = data.frame(
    "Sample_name" = colnames(count_data_centroid),
    "Literature" =  as.character(meta_data$Subtype)[cent_match],
    "Classification" =  as.character(groups),
    "Cor" = as.character(correlation),
    "P_value" = as.character(p_value),
    "Q_value" = as.character(q_value)
  )
  #classification_t[1:5,]
  
  write.table(
    classification_t,
    paste0(sep = "", collapse = "", c("~/MAPTor_NET/Results/Classification//Classification_IT",iteration,".tsv") ),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  to_be_excluded = as.character(classification_t$Sample_name[ as.double(as.character(classification_t$Q_value)) > 0.01]) # it0: 9k it1: 6526 # it2: 5287 # it3: 4536 it4: 4252 it5: 4148 it6: 4120 it7: 4116 it8: 4114
  
  ground_truth_vec <<- as.vector(classification_t$Classification[ !( as.character(classification_t$Sample_name) %in% to_be_excluded ) ])
  names(ground_truth_vec) = classification_t$Sample_name[ !( as.character(classification_t$Sample_name) %in% to_be_excluded ) ]
  nr_samples_new <<- length(ground_truth_vec)
  
  count_data = count_data[,!( colnames(count_data) %in% to_be_excluded)]
  source("~/MAPTor_NET//Scripts/Update_meta_info_table.R")
  
  iteration <<- iteration + 1
}
# 8547 7848 5278 5040

final_centroid_list = read.table(
  paste0(sep = "", collapse = "", c("~/MAPTor_NET//Results/Classification/Centroids_IT",iteration-1,".tsv") ),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

avg_importance = rowSums( abs((matrix(as.double(unlist(final_centroid_list[,3:ncol(final_centroid_list)])), ncol = ncol(final_centroid_list)-2)) ))
final_centroid_list = cbind( final_centroid_list, avg_importance)
final_centroid_list = final_centroid_list[order(final_centroid_list[,ncol(final_centroid_list)],decreasing = T),]

write.table(final_centroid_list,
            paste0(sep = "", collapse = "", c("~/MAPTor_NET//Results/Classification/Centroids_IT",iteration,".tsv") ),
            sep="\t",quote =F,row.names=F
)

meta_match = match(as.character(names(ground_truth_vec)), meta_info$Name, nomatch = 0)
meta_info$Classification[meta_match] = as.character(ground_truth_vec[ meta_match != 0])
meta_info$Classification = as.character(meta_info$Classification)
meta_info$Classification[ is.na(meta_info$Classification)  ] = ""
meta_info$Classification[ meta_info$Classification == ""  ] = "Not_sig"
write.table(meta_info, "~/MAPTor_NET//Misc/Meta_information.tsv",sep ="\t", row.names =T, quote =F)
