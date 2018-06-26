library("stringr")
library("ggplot2")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv", sep ="\t", stringsAsFactors = F, header = T)
meta_info = subset( meta_info, Included == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
colnames(meta_info) = str_replace_all(colnames(meta_info), pattern = "\\.AF8\\.","_")
f_mat = meta_info[ ,c("Subtype","Geschlecht","OS","Loc_Primaertumor","Treatment","Chemozyklen","Best_response","Vortherapie","Raucher","Alkohol","Anzahl_py")]
f_mat$OS = str_replace(f_mat$OS, pattern = ",",".")
f_mat$OS = as.double(f_mat$OS)
#f_mat$PFS = str_replace(f_mat$PFS, pattern = ",",".")
#f_mat$PFS = as.double(f_mat$PFS)

lin_mod = glm(
  data = f_mat,
  formula = OS ~ .#,
  #family=binomial(link='logit')
)
summary(lin_mod)


vis_mat = coef(lin_mod)
vis_mat = reshape2::melt(vis_mat)
vis_mat = cbind(rownames(vis_mat), vis_mat)
colnames(vis_mat) = c("Property","Impact")
vis_mat$Property = factor( vis_mat$Property, levels = vis_mat$Property[order(vis_mat$Impact)])
g = ggplot(
  vis_mat#,
  #aes( Property, Impact)
)
g = g + geom_bar(aes( x = Property, y = Impact),stat="identity")
g = g + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")
g

b0 <- lin_mod$coefficients[1]
b1 <- lin_mod$coefficients[2]

lin_vis_mat = data.frame(
    "Predicted" = fitted(lin_mod),
    "Reality" = f_mat$OS,
    "Subtype" = f_mat$Subtype,
    ""
)
cor(lin_vis_mat$Predicted, lin_vis_mat$Reality)
#lin_vis_mat$Predicted = scale(lin_vis_mat$Predicted)
#lin_vis_mat$Reality = scale(lin_vis_mat$Reality)

g_bench = ggplot(
  data = lin_vis_mat,
  aes( y =  Predicted,x = Reality))
g_bench = g_bench + geom_point( aes( color = lin_vis_mat$Subtype, size = 4))
g_bench = g_bench + geom_smooth(method = "lm")
#g_bench + xlab("Log2 variants per CCL-profile")
g_bench

### real expressen

r_mat = as.data.frame(matrix( as.double(as.character(unlist(pure_data) )), nrow = nrow(pure_data)))
colnames(r_mat) = colnames(pure_data)
rownames(r_mat) = rownames(pure_data)

r_mat = t(r_mat)
s_match = match(rownames(r_mat), meta_data$Name, nomatch = 0)
r_mat = r_mat[ s_match != 0,]

s_match = match(rownames(r_mat), meta_data$Name, nomatch = 0)
os_vec = meta_data$OS[s_match]
os_vec = str_replace(os_vec, pattern = ",",".")
os_vec = as.double(os_vec)
os_vec = os_vec[ s_match != 0]

r_mat = cbind(os_vec,r_mat)
colnames(r_mat)[1] = "OS"
rr_mat = as.data.frame(matrix( as.double( unlist(r_mat) ), nrow = nrow(r_mat), ncol = ncol(r_mat)))
colnames(rr_mat) = colnames(r_mat)
rownames(rr_mat) = rownames(r_mat)

rr_mat[1:5,1:5]
colnames(r_mat) = str_replace_all( colnames(r_mat), "\\.", "_" )

lin_mod = glm(
  data = data.frame(r_mat),
  formula = os_vec ~ . 
  #family=binomial(link='logit')
)
summary(lin_mod)
coff = coefficients(lin_mod)
coff = coff[!is.na(coff)]
coff_vis_mat = reshape2::melt(coff)
coff_vis_mat 

coff_vis_mat = cbind(rownames(coff_vis_mat), coff_vis_mat)
colnames(coff_vis_mat) = c("Property","Impact")
coff_vis_mat$Property = factor( coff_vis_mat$Property, levels = coff_vis_mat$Property[order(coff_vis_mat$Impact)])

coff_vis_mat = coff_vis_mat[abs(coff_vis_mat$Impact)> 50,]

g_bench = ggplot(
  data = coff_vis_mat,
  aes( y =  Property,x = Impact))
g_bench = g_bench + geom_bar(aes( x = Property, y = Impact),stat="identity")
g_bench + theme( axis.text.x = element_text(angle=45, vjust = .5))
