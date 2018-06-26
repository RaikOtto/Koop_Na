meta_data = data.frame("Subtype" = groups)
rownames(meta_data) = colnames(expr)

t = read.table("~/Koop_Klinghammer/Misc/Drugs.tsv",sep="\t",header =T, stringsAsFactors = F )
map = match(colnames(expr),t$Tumor.ID,nomatch = 0)

meta_data$Cetuximab[ map != 0] = t$Cetuximab[map]
