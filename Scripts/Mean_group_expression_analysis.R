# run group_allocation before!

cl_exp = result[,groups == "CL"]
ba_exp = result[,groups == "BA"]
ms_exp = result[,which( groups == "MS" )]

all_rowMeans = rowMeans( result )

cl_row_means = rowMeans(cl_exp)
ba_row_means = rowMeans(ba_exp)
ms_row_means = rowMeans(ms_exp) 

compare_data = c(
  cl_row_means,
  ba_row_means,
  ms_row_means
)
compare_data = data.frame(
  "data" = compare_data,
  "group" = factor(c(
    rep("CL",length(cl_row_means)),
    rep("BA",length(ba_row_means)),
    rep("MS",length(ms_row_means))
  ))
)
colnames(compare_data)

bartlett.test( x = compare_data$data, compare_data$group)
fligner.test(  x = compare_data$data, compare_data$group)

fit = lm(formula = compare_data$data ~ compare_data$group )
anova(fit)

t.test(ms_row_means,ba_row_means)

res_matrix = matrix(rep(0,9),nrow = 3)

colnames(res_matrix) = c("CL","BA","MS")
rownames(res_matrix) = c("CL","BA","MS")
res_matrix[1,1] = 1.0
res_matrix[2,2] = 1.0
res_matrix[3,3] = 1.0
res_matrix[1,2] = round(t.test(cl_row_means, ba_row_means)$p.value,3)
res_matrix[2,1] = round(t.test(cl_row_means, ba_row_means)$p.value,3)
res_matrix[1,3] = round(t.test(cl_row_means, ms_row_means)$p.value,3)
res_matrix[3,1] = round(t.test(cl_row_means, ms_row_means)$p.value,3)
res_matrix[2,3] = round(t.test(ba_row_means, ms_row_means)$p.value,3)
res_matrix[3,2] = round(t.test(ba_row_means, ms_row_means)$p.value,3)

### optical gimics


compare_data = cbind(
  round( cl_row_means, 3 ),
  round( ba_row_means, 3 ),
 round( ms_row_means, 3 )
)
compare_data = data.frame(
  "gene" = rownames(result),
  "data" = compare_data,
  "dif_cl_ba" = round( cl_row_means - ba_row_means, 3 ),
  "dif_cl_ms" = round( cl_row_means - ms_row_means, 3 ),
  "dif_ba_ms" = round( ba_row_means - ms_row_means, 3 )
)
colnames(compare_data) = c("HGNC","classic","basal","mesenchymal","dif_cl_minus_ba","dif_cl_minus_ms","dif_ba_minus_ms")
xlsx::write.xlsx2(compare_data,"~/Koop_Klinghammer/results/Averaged_gene_expression.xlsx")
write.table(compare_data, "~/Koop_Klinghammer/results/Averaged_gene_expression.tab",sep="\t",row.names = F, quote = F)
