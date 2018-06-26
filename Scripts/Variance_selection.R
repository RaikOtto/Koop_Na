
hgnc_list = hgnc_genes
hgnc_list_uni = unique(hgnc_list[!is.na(hgnc_list)])
len = as.double(length(hgnc_list_uni))
progress <<- 0.0
new_progress <<- progress

variance_selection = function( elem, input_matrix ){
  
  nr = as.double(which( elem == hgnc_list_uni ))
  
  new_progress <<- round( nr / len*100,0)
  if (new_progress != progress){
    print(  as.character( new_progress ) )
    progress <<- new_progress
  }
  
  mapping = as.character( hgnc_list ) == as.character( elem )
  mapping[ is.na(mapping)] = FALSE
  exprs_vals = matrix(
    as.double( 
      unlist( input_matrix[ mapping, ] )
    ),
    ncol = ncol(input_matrix)
  )
  exprs_var = apply( exprs_vals, MARGIN = 1, FUN = var  )
  selector  = which( exprs_var == max( exprs_var )  )[1]
  
  max_var_gene = matrix(
    as.double( exprs_vals[ selector, ] ),
    ncol = dim(input_matrix)[2]
  )
  colnames(max_var_gene) = colnames(input_matrix)
  
  return(max_var_gene)
}

max_list <<- c()
for (gene in hgnc_list_uni){
  var_match = which( hgnc_list == gene )
  row_var_max = which.max( apply(as.matrix(pure_data[ var_match,]),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1) )
  max_list = c( max_list, var_match[row_var_max])
}

res = pure_data[max_list,]

res = matrix(as.double(as.character(unlist(res))),ncol = ncol(pure_data) )
colnames(res) = colnames(pure_data)
rownames(res) = hgnc_list_uni[!is.na(hgnc_list_uni)]
expr = res
dim(expr)
