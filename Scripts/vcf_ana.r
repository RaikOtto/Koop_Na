library("stringr")
library("simpleaffy")
library("affy")
library("affyPLM")
library("RColorBrewer")
library("genefilter")
library("limma")
library("gplots")
library("arrayQualityMetrics")
library(GenomicFeatures)
library(VariantAnnotation)

setwd("~/Dropbox/PhD/Koop_Klinghammer/upload/VCF/")

# pre processing

subtypes  = read.table( "~/Dropbox/PhD/Koop_Klinghammer//subtype_classification.tab", sep ="\t", header = T  )
drug_data = read.table( "~/Dropbox/PhD/Koop_Klinghammer//drug_response_samples.tab" , sep = "\t", header = T )
mapping   = match(subtypes$Sample_ID , drug_data$Tumor.ID, nomatch = 0)
mapping_exists = mapping != 0

# start analysis

t = data.frame(
  "Cetuximab" =  drug_data$Cetuximab[mapping],
  "Subtype" = subtypes$Group[ mapping_exists ],
  stringsAsFactors = T
)

CL = drug_data$Cetuximab[ match( subtypes$Sample_ID[ subtypes$Group == "CL" ] , drug_data$Tumor.ID, nomatch = 0)  ]
MS = drug_data$Cetuximab[ match( subtypes$Sample_ID[ subtypes$Group == "MS" ] , drug_data$Tumor.ID, nomatch = 0)  ]
BA = drug_data$Cetuximab[ match( subtypes$Sample_ID[ subtypes$Group == "BA" ] , drug_data$Tumor.ID, nomatch = 0)  ]

res = c(CL,MS,BA)

# 

t2 = read.table("~/Dropbox/PhD/Koop_Klinghammer//drug_response_samples.tab",sep="\t",header=T)
t3 = as.matrix(t2)

heatmap.2(d,na.color="blue",dendrogram = "none", col = colorpanel(100,"green","white","red"),trace="none")

setwd("~/Dropbox/PhD/Koop_Klinghammer//upload/Cel_files/")

celFiles = list.celfiles( "." )

rawdata = ReadAffy( filenames= celFiles )
eset_ori = rma( rawdata )
eset = eset_ori

#qc_res = qc(rawdata_k)
#degrade<-AffyRNAdeg(rawdata_k)
#plotAffyRNAdeg(degrade)

annot_data    = read.table("../../Generic_Biomarker_mRNA_Pipeline/HG-U133_Plus_2.na34.annot.csv",sep=",",header=T)

mapping       = match( rownames( eset ), annot_data$Probe.Set.ID )
hgnc_genes    = as.character( annot_data$Gene.Symbol[ mapping ] )
ensembl_genes = as.character( annot_data$Ensembl[ mapping ] )
entrez_genes  = as.character( annot_data$Entrez.Gene[ mapping ] )
pathway       = as.character( annot_data$Pathway[ mapping ] )
omim          = as.character( annot_data$OMIM[ mapping ] )

pData( eset ) = cbind( pData( eset ), rownames( pData( eset ) ) )
#colnames ( pData(eset_k) ) = c("index","sample")

drug_t = read.table("~/Dropbox/PhD/Klinghammer/drug_sens.tab",sep="\t", header=T)
s=c("9619","9876","9897","10110","10114","10159","10309","10321","10379","10511","10621","10632","10847","10913","10924","10927","10960","11097","11218","10980A","10980B","11269A","11269B","11437A","11437B","11452","11482A","11142","11143","10883","11204A","11857B","11841","11204B","11527A")
rownames(drug_t) = s

f = function( vec ){ return( median( sort(vec[ !is.na(vec) ]) ) ) }

ord_pat= order(apply( drug_t, FUN= f, MARGIN=1 ))

boxplot(t(drug_t)[,ord_pat],xlab="study_id",xaxt="n", main="Drug sensitivity per study_id",ylab="Tumor growth")

axis(1, at= 1:35, labels = FALSE)
text(y=-1,x = 1:35, labels = s[ord_pat], srt = 45, pos=1,xpd=T)
#axis(1,at=1:35, labels = ord_pat,srt=45)

ord_drug = order(apply( drug_t, FUN= f, MARGIN=2 ))

boxplot( drug_t[,ord_drug],xlab="drug",xaxt="n", main="Drug sensitivity per drug",ylab="Tumor growth")
axis(1, at= 1:8, labels = FALSE)
text(y=-1,x = 1:8, labels = colnames(drug_t)[ ord_drug ], srt = 45, pos=1,xpd=T)
#axis(1,at=1:8, labels = ord_drug,srt=45)

# eset <-threestep(  rawdata_k, background.method="RMA.2", normalize.method="quantile", summary.method="median.polish")

m2 = match( str_replace(colnames(eset),".CEL",""), rownames( drug_t ) )

# 29940

hgnc_genes[which( entrez_genes == "29940" ) ]
m = which( entrez_genes == "29940" ) # SART2

d = cbind( exprs(eset)[m[1],], drug_t$Docetaxel[m2])
d = d[ order(d[,2]  ),]

d2 = d [ -which( rownames(d) == "10321.CEL" ) ,]
plot( d2[,2], d2[,1] )

# 7798

hgnc_genes[which( entrez_genes == "7798" ) ]
m = which( entrez_genes == "7798" ) # LUP1

d = cbind( exprs(eset)[m[1],], drug_t$Docetaxel[m2])
d = d[ order(d[,2]  ),]
cor( d[,1], d[,2], method = "spearman", use = "pairwise.complete.obs"  )

plot( d[,2], d[,1] )

fit = lm ( d[,2] ~ d[,1]  )
abline( fit, col="red")

d2 = d [ -which( rownames(d) == "10321.CEL" ) ,]
cor( d2[,2], d2[,1], method = "spearman", use = "pairwise.complete.obs"  )
plot( d2[,2], d2[,1] )

fit = lm ( d2[,1] ~ d2[,2]  )
abline( fit, col="red")

# 5243

hgnc_genes[which( entrez_genes == "5243" ) ]
m = which( entrez_genes == "5243" ) # ABCB1

d = cbind( exprs(eset)[m[1],], drug_t$Docetaxel[m2])
d = d[ order(d[,2]  ),]
cor( d[,1], d[,2], method = "spearman", use = "pairwise.complete.obs"  )

plot( d[,2], d[,1] )

fit = lm ( d[,2] ~ d[,1]  )
abline( fit, col="red")

d2 = d [ -which( rownames(d) == "10321.CEL" ) ,]
cor( d2[,2], d2[,1], method = "spearman", use = "pairwise.complete.obs"  )
plot( d2[,2], d2[,1] )

fit = lm ( d2[,1] ~ d2[,2]  )
abline( fit, col="red")