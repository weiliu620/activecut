import sys
import os
import numpy as np
import nibabel as nib
import subprocess


def geteigenimage(eigenvec, eigenimage_prefix, coord, mask_file, n_comp):
    """
    Map the princiapl components to 3D brain space. eigenvec is a sample x D matrix. 
    """
    mask = nib.load(mask_file)
    (nvars, neigs) = eigenvec.shape

    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    for eig_idx in range(n_comp):
        (ncoord, dummy) = coord.shape
        D = np.zeros(mask.shape)
        D[:] = eigenvec[:,eig_idx].min()
        print 'eig vec {0}, min: {1:.5f}, max: {2:.5f}'.format(eig_idx, eigenvec[:,eig_idx].min(), eigenvec[:,eig_idx].max())
        for coord_idx in range(ncoord):
            vox_idx = np.dot(S, coord[coord_idx,:] - offset)
            vox_idx = tuple(np.around(vox_idx).astype(int))
            
            draw_roi(D, vox_idx, radius, eigenvec[coord_idx, eig_idx])

        # save this image.
        eigenimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
        eigenimage.set_data_dtype('float32')
        eigenimage.to_filename(eigenimage_prefix + str(eig_idx).zfill(3) + '.nii.gz')
        print '{0} saved. '.format(eigenimage_prefix + str(eig_idx).zfill(3) + '.nii.gz')

        
def draw_roi(D, vox_idx, radius, value):
    """
    draw a ball-shape ROI around the seed index
    """

    cand = [(i, j, k) for i in range(vox_idx[0] - radius, vox_idx[0] + radius + 1) for j in range(vox_idx[1] - radius, vox_idx[1] + radius + 1) for k in range(vox_idx[2] - radius, vox_idx[2] + radius + 1) if (i-vox_idx[0])**2 + (j-vox_idx[1])**2 + (k-vox_idx[2])**2 <= radius**2]
    for v in cand:
        D[v] = value
        
        
def draw_sample_pos(coord, mask_file, out_file):
    """
    draw sample location on a 3d template volume
    """
    mask = nib.load(mask_file)
    (ncoord, dummy) = coord.shape
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    D = np.zeros(mask.shape, dtype = 'uint8')
    for coord_idx in range(ncoord):
        vox_idx = np.dot(S, coord[coord_idx,:] - offset)
        vox_idx = tuple(np.around(vox_idx).astype(int))

        cand = [(i, j, k) for i in range(vox_idx[0] - radius, vox_idx[0] + radius + 1) for j in range(vox_idx[1] - radius, vox_idx[1] + radius + 1) for k in range(vox_idx[2] - radius, vox_idx[2] + radius + 1) if (i-vox_idx[0])**2 + (j-vox_idx[1])**2 + (k-vox_idx[2])**2 <= radius**2]
        for v in cand:
            D[v] = D[v] + 1
        
    I = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
    I.set_data_dtype('uint8')
    I.to_filename(out_file)

def draw_label_map(labels, out_file, coord, mask_file):
    """
    draw the clustering label map to the 3d volume.
    """
    mask = nib.load(mask_file)
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    (ncoord, dummy) = coord.shape
    D = np.zeros(mask.shape, dtype = 'uint8')

    for coord_idx in range(ncoord):
        vox_idx = np.dot(S, coord[coord_idx,:] - offset)
        vox_idx = tuple(np.around(vox_idx).astype(int))

        cand = [(i, j, k) for i in range(vox_idx[0] - radius, vox_idx[0] + radius + 1) for j in range(vox_idx[1] - radius, vox_idx[1] + radius + 1) for k in range(vox_idx[2] - radius, vox_idx[2] + radius + 1) if (i-vox_idx[0])**2 + (j-vox_idx[1])**2 + (k-vox_idx[2])**2 <= radius**2]
        for v in cand:
            D[v] = labels[coord_idx] + 1
            
    # save this image.
    outimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
    outimage.set_data_dtype('uint8')
    outimage.to_filename(out_file)
    print '{0} saved. '.format(out_file)        
    
def exp_hist_map(pacall, coord_file, mask_file, out_file):
    """
    draw the number of expressed genes at each regions. Normalized to 1.
    """
    mask = nib.load(mask_file)
    rsum = pacall.sum(axis = 0).astype(float)
    coord = np.loadtxt(coord_file)
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    
    (n_genes, n_samples) = pacall.shape
    D = np.zeros(mask.shape)
    for sample_idx in range(n_samples):
        vox_idx = np.dot(S, coord[sample_idx,:] - offset)
        vox_idx = tuple(np.around(vox_idx).astype(int))
        draw_roi(D, vox_idx, radius, rsum[sample_idx] / n_genes)

    # save this image.
    eigenimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
    eigenimage.set_data_dtype('float32')
    eigenimage.to_filename(out_file)
    print '{0} saved. '.format(out_file)


def map_to_brain(vec, coord_file, mask_file, out_file):
    """
    A general routine to map a vector to the ROIs in the brain.
    """

    print 'vector min: {0}. max = {1}'.format(vec.min(), vec.max())
    coord = np.loadtxt(coord_file)
    mask = nib.load(mask_file)
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    
    (ncoord, dummy) = coord.shape
    D = np.zeros(mask.shape)
    for coord_idx in range(ncoord):
        vox_idx = np.dot(S, coord[coord_idx,:] - offset)
        vox_idx = tuple(np.around(vox_idx).astype(int))
        draw_roi(D, vox_idx, radius, vec[coord_idx])

    # save this image.
    eigenimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
    eigenimage.set_data_dtype('float32')
    eigenimage.to_filename(out_file)
    print '{0} saved. '.format(out_file)

def map_to_brain_color(vec, coord, mask_file, out_file):
    """
    map Px3 array to RGB brain map.
    """
    mask = nib.load(mask_file)
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)
    (n_tisamples, dummy) = coord.shape
    
    for eig_idx in range(0,3):
        # convert to [0, 255] for rgb visualizattion.
        inten_min = np.percentile(vec[:, eig_idx], 0.1)
        inten_max = np.percentile(vec[:, eig_idx], 99.9)
        # outliers are clipped.
        vec[:,eig_idx] = 255 * (vec[:,eig_idx] - inten_min) / (inten_max - inten_min)
        vec[vec[:,eig_idx] > 255, eig_idx] = 255
        vec[vec[:,eig_idx] < 0,   eig_idx] = 0
        D = np.zeros(mask.shape)
        for tisample_idx in range(n_tisamples):
            vox_idx = np.dot(S, coord[tisample_idx,:] - offset)
            vox_idx = tuple(np.around(vox_idx).astype(int))
            draw_roi(D, vox_idx, radius, vec[tisample_idx, eig_idx])
            # print 'eig_idx: {0}, tisample_idx: {1}, intensity: {2}'.format(eig_idx, tisample_idx, vec[tisample_idx, eig_idx])
            
        # save this image.
        eigenimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
        eigenimage.set_data_dtype('float32')
        eigenimage.to_filename('channel_' + str(eig_idx).zfill(3) + '.nii.gz')
        print '{0} saved.'.format('channel_' + str(eig_idx).zfill(3) + '.nii.gz')

        # convert to nrrd.
        subprocess.call(['convertITKformats', 'channel_' + str(eig_idx).zfill(3) + '.nii.gz', 'channel_' + str(eig_idx).zfill(3) + '.nhdr'])
        
    subprocess.call(['unu', 'join', '-i', 'channel_000.nhdr', 'channel_001.nhdr', 'channel_002.nhdr', '-a', '0', '-incr', '-o', out_file ] )
    
    os.unlink('channel_' + str(0).zfill(3) + '.nii.gz')
    os.unlink('channel_' + str(1).zfill(3) + '.nii.gz')
    os.unlink('channel_' + str(2).zfill(3) + '.nii.gz')
    os.unlink('channel_' + str(0).zfill(3) + '.nhdr')
    os.unlink('channel_' + str(1).zfill(3) + '.nhdr')
    os.unlink('channel_' + str(2).zfill(3) + '.nhdr')
    os.unlink('channel_' + str(0).zfill(3) + '.raw.gz')
    os.unlink('channel_' + str(1).zfill(3) + '.raw.gz')
    os.unlink('channel_' + str(2).zfill(3) + '.raw.gz')
