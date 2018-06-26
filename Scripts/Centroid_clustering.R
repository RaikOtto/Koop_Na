
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
      pv[ j ] = corj$p.value
    }
    
    maxk = which.max(c.cor)
    pub_cor <<- rbind(pub_cor,c.cor)
    group = names( maxk )
    vcor = rbind( vcor, c( colnames( vd )[ i ], group, c.cor[ maxk ], pv[ maxk ] ) )
    
    if( pv[ maxk ] < .05 ){
      vclass[ colnames( vd )[ i ] ] = group
    }
  }
  
  
  return( list( overlap.gene = gene.sig, cluster = vclass, correlation = vcor ) )
}

#pub_cor <<- matrix( as.double(), ncol = length( c("Cl1","Cl2_5_6","Cl3_4")  ) )

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )
expr2bc = centroid2expr( balanced.centroid, result_mat )
corrected = p.adjust( expr2bc$correlation[,4], "BH" )
expr2bc$correlation = cbind( expr2bc$correlation, corrected  )

shiny_mat = apply( expr2bc$correlation[,c(3,4,5)], FUN = as.double, MARGIN = 2  )
shiny_mat[, 1] = round( shiny_mat[,1], 2)
expr2bc$correlation[,c(3:5)] = shiny_mat
colnames(expr2bc$correlation) = c("Sample","Subtype","Correlation","P_Value","Q_value")

write.table( 
  expr2bc$correlation,
  file = "~/Koop_Klinghammer//Sub_typ_class_new_01_11.tab",
  quote = F,
  row.names = F,
  sep = "\t"
)
