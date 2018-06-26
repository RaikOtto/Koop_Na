library("stringr")
library("caret")

r_mat = read.table("~/Koop_Klinghammer/Data/Pure_data.05_06_2018.tsv",sep="\t", header = T)
colnames(r_mat) = str_replace_all( colnames(r_mat), pattern = "^X", "")
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv", sep ="\t", stringsAsFactors = F, header = T)
meta_info = subset( meta_info, Included == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info = subset( meta_info, Sig == TRUE)
meta_info$OS = as.double( str_replace_all(meta_info$OS, pattern = ",", "\\."))

###

s_match = match(rownames(t(r_mat)), meta_info$Name, nomatch = 0)
r_mat = r_mat[s_match != 0,]
meta_data = meta_info[ s_match, ]

f_mat = meta_info[ ,c("Subtype","Geschlecht","OS","Loc_Primaertumor","Treatment","Chemozyklen","Best_response","Vortherapie","Raucher","Alkohol","Anzahl_py")]
f_mat$OS = str_replace(f_mat$OS, pattern = ",",".")
f_mat$OS = as.double(f_mat$OS)


r_mat = cbind( as.character( rownames(r_mat) ) , r_mat)
colnames(r_mat)[1] = "Predictors"

r_mat = rbind( as.character( c( "", meta_data$Subtype ) ), r_mat)
rownames(r_mat)[1] = "Subtype"

r_mat = t(r_mat)

#lambda <- 10^seq(10, -2, length = 100)

splitIndex = caret::createDataPartition( 1:nrow(r_mat), p = .75, list = FALSE, times = 1)
trainDF = r_mat[ splitIndex,]
testDF  = r_mat[-splitIndex,]

predictorsNames = colnames(r_mat)[2:(ncol(r_mat))]

tune_grid = tuneGrid = expand.grid(
  .alpha=1,
  .lambda=seq(0, 100, by = 0.1)
)

objControl <- caret::trainControl(method='cv', number=3, returnResamp='none', classProbs = TRUE)

trainControl = caret::trainControl(
  method = "cv",
  number = 10,
  summaryFunction = prSummary,
  classProbs = T
)

prop.table(table(r_mat[,1]))

modelFit = caret::train(
  trainDF[,predictorsNames],
  trainDF[,"Subtype"],
  #pure_data,
  #r_mat[,"Subtype"],
  method='svmLinearWeights',
  trControl=objControl,  
  metric = "ROC"#,
  #preProc = c("center", "scale")
)

predictions <- predict(object = modelFit, testDF[,predictorsNames], type='raw')
table(as.character(predictions), testDF[,1])
#auc <- pROC::roc( testDF[,outcomeName], predictions)
#print(auc$auc)

d=varImp(modelFit,scale=T)
imp_vars = d$importance
vis_vars = imp_vars[order(imp_vars[,1], decreasing = T),]
names(vis_vars) = rownames(imp_vars)[order(imp_vars[,1], decreasing = T)]

plot( vis_vars[1:10] )

aggregate( pure_data[rownames(pure_data) == "SPRR3",], by = list(meta_data$Subtype), FUN = mean)
###

#create test and training sets

vimp <- varImp(modelFit, scale=F)
results <- data.frame(row.names(vimp$importance),vimp$importance$Overall)
results$VariableName <- rownames(vimp)
colnames(results) <- c('VariableName','Weight')
results <- results[order(results$Weight),]
results <- results[(results$Weight != 0),]

par(mar=c(5,15,4,2)) # increase y-axis margin. 
xx <- barplot(results$Weight, width = 0.85, 
              main = paste("Variable Importance -",outcomeName), horiz = T, 
              xlab = "< (-) importance >  < neutral >  < importance (+) >", axes = FALSE, 
              col = ifelse((results$Weight > 0), 'blue', 'red')) 
axis(2, at=xx, labels=results$VariableName, tick=FALSE, las=2, line=-0.3, cex.axis=0.6)  

###

write.table(r_mat, "~/Koop_Klinghammer/Results/Linear_regression/r_mat.tsv",quote = F, row.names = F, sep ="\t")

### SVM

library(e1071)
data(iris)
attach(iris)

train_data = data.frame( meta_data$Subtype, t(r_mat) )
train_data[1:5,1:5]
colnames(train_data)[1] = "Subtype"

svm_tune <- tune(svm, Subtype ~ . , data = train_data,
                 ranges = list(epsilon = seq(0,1,0.01), cost = 2^(2:9))
)
print(svm_tune)

model <- svm( Subtype ~ . ,train_data  , probability = TRUE, kernel = "radial" , nu = 1)
pred_prob <- predict(model, train_data, decision.values = TRUE, probability = TRUE)
as.character(pred_prob[1:ncol(pure_data)])

meta_data$Subtype = as.character(pred_prob[1:ncol(pure_data)])

min_vec = apply(matrix(( as.double(as.character(unlist(pred_prob[3])))), ncol = 3), FUN = function(vec){
    return( which.min(vec) )
  },
  MARGIN = 1)
