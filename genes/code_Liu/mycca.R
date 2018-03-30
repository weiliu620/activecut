library("CCA")
library("R.matlab")
X = readMat("/usr/sci/projects/genes/human/allsubjects/expressions.mat")
Y = readMat("/usr/sci/projects/genes/human/allsubjects/fmri/NY_TRT/session1_allsubjects.mat")
rcc(X$exp.allsub, Y$fmri, 0.1, 0.2)

