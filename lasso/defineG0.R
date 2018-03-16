defineG0 <- function(linear2sub, eps)
  ## This function gives a G0 prior network given voxels'
  ## linear-to-subscript projection matrix, and epsilon distance.
### linear2sub: Px3 matrix, where P is number of nodes (voxels). Each line is node p's voxel coordinates.
### eps: a distance. Within eps, two voxel (x1, y1, z1) and (x2, y2, z2) will be defined as neighbors, and G(voxel1, voxel2) will be 1.
### Output: a PxP matrix, containing either 1 or zero.

  {
    nNodes = nrow(linear2sub) #number of nodes.
    G0 = matrix(0, nrow = nNodes, ncol = nNodes)

    for (i in 1:nrow(G0) ) {
      for (j in 1:ncol(G0) ) { 
        ## node i and j's distance is bigger than eps?
        isub = linear2sub[i,] # i's subscript
        jsub = linear2sub[j,]
        if ( (isub[1] - jsub[1])^2 +  (isub[2] - jsub[2])^2 +  (isub[3] - jsub[3])^2 < eps^2) {
          G0[i,j] = 1
        }
        else {
          G0[i,j] = 0
        } 
      } # for j
      
    } # for i

    diag(G0) = 0
    G0 
  }
