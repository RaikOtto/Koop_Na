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