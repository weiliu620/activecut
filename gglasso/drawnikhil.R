betan = read.table("gridbetaactivesetFile", fill=TRUE, col.names=seq(1,20))
dev.new()
ts.plot(betan, main="Nikhil's solution")



# Hastie's solution
library("glmnet")
x = as.matrix(read.table("xx.txt"))
y = as.matrix(read.table("yy.txt"))
allLambda = as.matrix(read.table("alllambda.txt"))
fit = glmnet(x, y, lambda=allLambda, standardize=FALSE, type.gaussian="covariance")
dev.new()
plot(fit)  



