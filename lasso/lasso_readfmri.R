lasso.readfmri <- function(imageDir,
                           maskFile,
                           roiAtlas,
                           maskThresh,
                           beginSliceIdx,
                           endSliceIdx)
  ### imageDir: dir of source fmri images.
  ### maskFile: filename of brain mask images. Will be with absolute/relative path.
  ### beginSliceIdx and endSliceIdx: to extract a few slices since volume is too big.

  ### The function returns a list:
  ### Y: original NPxT data matrix, where N is # of sub, P is # of nodes(voxels), and T is time series length.
  ### sub2linear: map subscripts (x, y, z) to linear index p. (row in the PxT matirx)
  ### linear2sub: map p to (x, y, z)
  ### filenames: all file names of the fMRI files.
  ### beginSliceIdx and endSliceIdx: extracted slices.

  ## The region of interest definition (values are from fsl atlas Harvard cortex 2mm)
  ## PCC: Posterior Cingular Cortex. value = 30.
  ## MTL: Medial Temporal lobes. (approximated). Value = 21 or 20
  
  
  {
    # load libraries    
    library("Rniftilib")

    allfiles <- list.files(full.names = TRUE, path = imageDir)


    ## We need mask file, so only voxels in the brain are
    ## extracted. Although all subjects image are registered to same
    ## MNI space, they have different segmentation resutlts and have
    ## their own mask gray matter mask file. It may be good to have a
    ## single mask file that can be used for all images, but there is
    ## no good method to quicly do that  now.

    ## So for this code, I use subject 1's mask for all the subjects. This is,
    ## for sure, not the right way, but should not have impact on following test.

    maskStruct <- nifti.image.read(maskFile)
    roiAtlasStruct <- nifti.image.read(roiAtlas)

    ## Define a M by T matrix Y, where M = NP. P is number of nodes
    ## (voxels), and N is number of subjects. T is time series length.

    ## But we don't know P now, so have to define an empty matrix.
    fmriStruct <- nifti.image.read(allfiles[1])
    Y <- matrix(0, 0, fmriStruct$dim[4])

    # mask out points where any subject has zero time series.
    for (fileIdx in 1:length(allfiles)) {
      cat("file ", fileIdx, ": ", allfiles[fileIdx], "\n")
      
      for (x in 1:maskStruct$dim[1]) {
        for (y in 1:maskStruct$dim[2]) {
          for (z in beginSliceIdx:endSliceIdx) {
            if (sum(fmriStruct[x, y, z,]) == 0) {
              maskStruct[x, y, z] = 0
            }

          } # z
        } # y
      } # x
    }
    
    ## save modified mask image.
    ## nifti.set.filenames(maskStruct, "modifiedmask.nii")
    ## nifti.image.write(maskStruct)


      
### Compute mapping matrix sub2linear and linear2sub.
      
    ## each coordinate (x, y, z) of sub2linear has a value is the
    ## linear index.
    sub2linear <- array(0, dim = maskStruct$dim)
      
    ## we don't know how many gray matter voxels, so init as an empty
    ## matrix. n'th row of the matrix is 3-element vector (x, y, z),
    ## which is the (voxel coordinates) of voxel (or nodte) n.
    linear2sub <- matrix(0, 0, 3)

    ## make en empty matrix to save the region of interest
    ## (ROIs). Will be PxK, where P is # of nodes and K is # of
    ## ROIs. Not I just extract Posterior Cingular Cortex (PCC) so
    ## there is only one column.
    roiMask <- matrix(0, 0, 2)
    colnames(roiMask) <- c("PCC", "MTL")
      
    linearIdx = 0;
      for (x in 1:maskStruct$dim[1]) {
        for (y in 1:maskStruct$dim[2]) {
          for (z in beginSliceIdx:endSliceIdx) {
            if (maskStruct[x, y, z] > maskThresh) {
              ##
              linearIdx <- linearIdx + 1
              sub2linear[x, y, z] <- linearIdx

              ## Save current coordinates.
              linear2sub <- rbind(linear2sub, c(x, y, z))

              ## build roi mask.
              ## map to world coordinates.

              worldCoord <- fmriStruct$qto.xyz %*% as.matrix(c(x, y, z, 1))
              ## world coordinates map to voxel coordinates in std
              ## space.
              roiAtlasCoord <- roiAtlasStruct$qto.ijk %*% worldCoord

              isPCC <- roiAtlasStruct[roiAtlasCoord[1], roiAtlasCoord[2], roiAtlasCoord[3]] == 30
              isMTL <- roiAtlasStruct[roiAtlasCoord[1], roiAtlasCoord[2], roiAtlasCoord[3]] == 20 || roiAtlasStruct[roiAtlasCoord[1], roiAtlasCoord[2], roiAtlasCoord[3]] == 21
              roiMask <- rbind(roiMask, c(isPCC, isMTL ) )
            }
            
          } # z
        } # y
      } # x

    ## Extract each time series and build a P by T matrix Y.p. P is
    ## number of voxels in the brain's gray matter mask, and is the
    ## length of time series. N is the subject number.

    for (fileIdx in 1:length(allfiles)) {
      print(allfiles[fileIdx])
      fmriStruct <- nifti.image.read(allfiles[fileIdx])

      ## Define an empty 0 by T matrix, we'll append new data at end.
      Y.p = matrix(0, 0, fmriStruct$dim[4]);
      for (x in 1:maskStruct$dim[1]) {
        for (y in 1:maskStruct$dim[2]) {
          for (z in beginSliceIdx:endSliceIdx) {
            
            if (maskStruct[x, y, z] > maskThresh) {
              Y.p <- rbind(Y.p, fmriStruct[x, y, z,])

            } # end of if
          } # for z
        } # for y
      } # for x

      # add row names
      myrowname <- paste("sub", as.character(fileIdx), "node", 1:dim(Y.p)[1])
      rownames(Y.p) <- myrowname

      # Now I got all data for current subject. Concatenate to Y.
      Y <- rbind(Y, Y.p)
    }

    ## Do PCA on the NP by T matrix Y. NP is number of replicates and
    ## T is number of dimensions. So this fit the requirement of prcomp func.

    pcaObj <- prcomp(Y)
    
    return (list(Y = Y,
                 pcaObj = pcaObj,
                 sub2linear = sub2linear,
                 linear2sub = linear2sub,
                 roiMask = roiMask,
                 filenames = allfiles,
                 beginSliceIdx = beginSliceIdx,
                 endSliceIdx = endSliceIdx,
                 nSub = length(allfiles),
                 nNodes = nrow(linear2sub),
                 tsLength = fmriStruct$dim[4]))

}
