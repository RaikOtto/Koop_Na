library("stringr")
library("limma")

design <- model.matrix(~0 + meta_data$relapse_dichotomic)
colnames(design) = c("Healthy","No", "Yes")

vfit <- lmFit(pure_data,design)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#plotSA(efit)

contr.matrix = makeContrasts( contrast = Yes - No ,  levels = design )
vfit <- contrasts.fit( vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))

result_t = topTable( efit, coef = "contrast", number  = nrow(pure_data), adjust  ="none", p.value = 1, lfc = 0)
result_t$hgnc_symbol = rownames(result_t)
colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
result_t = result_t[order(result_t$P_value, decreasing = F),]
result_t$Log_FC = round(result_t$Log_FC, 1)
result_t$Average_Expr = round(result_t$Average_Expr, 1)
result_t = result_t[order(result_t$Log_FC,decreasing = T),]

write.table("~/Koop_Na/Results/Prelim_ana_26_06_2018/Dif_exp_dichotomic.tsv", x = result_t, sep = "\t", quote = F, row.names = F)

