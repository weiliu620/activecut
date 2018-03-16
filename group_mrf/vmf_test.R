## use movMF package to test my synthetic data.
library("Rniftilib")
fmriFile = "/home/sci/weiliu/dataset/synthetic_fmri_for_vmf_R/sub01.nii"
initLabelFile = "/home/sci/weiliu/dataset/synthetic_fmri_for_vmf_R/best_4_clusters_syn.nii"
trueLabelFile = "/home/sci/weiliu/dataset/synthetic_fmri_for_vmf_R/truelabel_sub01.nii"

initStruct <- nifti.image.read(initLabelFile)
fmriStruct <- nifti.image.read(fmriFile)
truelabelStruct <- nifti.image.read(trueLabelFile)

## construct data matrix.
totalPts = initStruct$dim[1] * initStruct$dim[2] * initStruct$dim[3];

linear2sub <- matrix(0, totalPts, 3)
sub2linear <- array(0, dim = initStruct$dim)
fmriMat = matrix(0, totalPts, fmriStruct$dim[4])
initLabel = matrix(0, totalPts, 1)
trueLabel = matrix(0, totalPts, 1)

linearIdx = 0;
for (x in 1:initStruct$dim[1]) {
  for (y in 1:initStruct$dim[2]) {
    for (z in 1:initStruct$dim[3]) {
      if (initStruct[x, y, z] > 0) {
        linearIdx <- linearIdx + 1
        fmriMat[linearIdx,] = fmriStruct[x, y, z,]
        initLabel[linearIdx,] = initStruct[x, y, z]
        trueLabel[linearIdx,] = truelabelStruct[x, y, z]
        sub2linear[x, y, z] <- linearIdx
        linear2sub[linearIdx,] <- c(x, y, z)
        print(linearIdx)
      }
    }
  }
}

# truncate
fmriMat = fmriMat[1:linearIdx,]
initLabel = initLabel[1:linearIdx,]
linear2sub = linear2sub[1:linearIdx,]
trueLabel = trueLabel[1:linearIdx,]

print("done with collecting data.")
set.seed(123)
myvmmodel = movMF(fmriMat, 4, control=list(maxiter=20, start=initLabel, verbose=TRUE))
predictedLabel = predict(myvmmodel)

for (n in 1:linearIdx) {
  x = linear2sub[n,1]
  y = linear2sub[n,2]
  z = linear2sub[n,3]
  
  initStruct[x, y, z] = as.numeric(predictedLabel[n])
}

nifti.set.filenames(initStruct, "predictedLabelMap")
nifti.image.write(initStruct)




  

