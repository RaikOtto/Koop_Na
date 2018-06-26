require('NanoStringNorm')
library("stringr")

meta_info = read.table("~/Koop_Na//Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)

raw_data  = read.markup.RCC( rcc.path = "~/Koop_Na/Data/Raw_data/", rcc.pattern = "*.RCC")

### tmp
sample_names = names(raw_data$header)
sample_names = str_replace(sample_names, pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")

meta_match = match( sample_names, meta_info$Raw_Name, nomatch = 0)
colnames(raw_data$x)[-seq(3)] = meta_info$Name[meta_match]

### normalization

#norm.comp.results.test = norm.comp(raw_data, verbose = T)
eset = NanoStringNorm::NanoStringNorm( 
  raw_data,
  CodeCount.methods = "sum",
  Background.methods = "mean.2sd",
  SampleContent.methods = "housekeeping.sum",
  OtherNorm.methods = "vsn",
  take.log = T,
  round.values = T,
  return.matrix.of.endogenous.probes = F,
  #traits = meta_data,
  verbose = T
)

m    = matrix( as.character(unlist( eset$normalized.data)), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2])
info = m[,seq(3)]
data = matrix( as.double(m[,-seq(3)]), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2]-3)
data = round(data,1)

rownames(data) = rownames(eset$normalized.data)
col_labels = str_replace( colnames(eset$normalized.data)[-seq(3)], pattern = "^X", "") 

colnames(data) = col_labels
res = cbind( info,data )
res = cbind(rownames(eset$normalized.data), res)

pure_data = as.character(res)[-seq(dim(eset$normalized.data)[1]*4)]
pure_data = matrix( as.double( pure_data ), nrow = dim(eset$normalized.data)[1] )
rownames( pure_data ) = info[ ,2 ]
pure_data = pure_data[ info[,1] == "Endogenous"  ,]

s_match = match( col_labels, as.character( meta_info$Name), nomatch = 0)
meta_data = meta_info[s_match,]
rownames(meta_data) = meta_data$SampleID
colnames(pure_data) = meta_data$SampleID

### optional normalization

design <- model.matrix(~0 + meta_data$relapse_dichotomic)
colnames(design) = c("Healthy","No", "Yes")

DGE = edgeR::DGEList(pure_data)
DGE = edgeR::calcNormFactors(DGE,method =c("TMM"))
v = limma::voom(DGE,design)

#pure_data = DGE$counts
#pure_data = v$E
boxplot(pure_data)

#write.table(pure_data, "~/Koop_Klinghammer/Data/Pure_data.05_06_2018.tsv", quote= F, row.names = T, sep = "\t")
