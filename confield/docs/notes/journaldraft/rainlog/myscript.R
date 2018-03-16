kmeans = read.table('kmeans.log')
ncuts = read.table('ncuts.log')
groupmrf_1g = read.table('groupmrf_test1g.log')

s12 <- cbind(kmeans = kmeans[,2], ncuts = ncuts[,2], groupmrf = groupmrf_1g[,2])
s13 <- cbind(kmeans = kmeans[,3], ncuts = ncuts[,3], groupmrf = groupmrf_1g[,3])
s23 <- cbind(kmeans = kmeans[,4], ncuts = ncuts[,4], groupmrf = groupmrf_1g[,4])


t.test(s12[,"groupmrf"], s12[,"kmeans"], paired=TRUE)
t.test(s13[,"groupmrf"], s12[,"kmeans"], paired=TRUE)
t.test(s23[,"groupmrf"], s12[,"kmeans"], paired=TRUE)

dev.new(width=10, height=5)
par(mfcol=c(1,3), cex.axis=1.5,cex.main=1.5)
mylim=c(0.74, 0.9)
boxplot(s12, main = "Rand index between session 1 and 2", col="lightgray", ylim=mylim)
boxplot(s13, main = "Rand index between session 1 and 3", col="lightgray", ylim=mylim)
boxplot(s23, main = "Rand index between session 2 and 3", col="lightgray", ylim=mylim)
dev.copy2pdf(file="boxplot.pdf")

