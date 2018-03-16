library("grplasso")
data(splice)
contr <- rep(list("contr.sum"), ncol(splice)-1)
names(contr) <- names(splice)[-1]
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20, contrasts = contr, center = TRUE, standardize = TRUE)
