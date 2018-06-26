# unsupervized data analysis

par(mfrow=c(1,2))

library(RColorBrewer)

sample_names = colnames(result)
nsamples = length(sample_names)

col = brewer.pal(nsamples, "Paired")

plot(
  density(result[,1]),
  col=col[1],
  lwd=2,
  ylim=c(0,0.21),
  las=2, 
  main="",
  xlab=""
)
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:nsamples){
  den <- density(result[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample_names, text.col=col, bty="n")

limma::plotMA(result)

result_save = result

