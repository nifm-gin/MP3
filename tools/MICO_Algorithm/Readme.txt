This package includes the Matlab code that implements the method for tissue segmentation and bias field correction in Chunming Li et al's paper:
    "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation", Magnetic Resonance Imaging, vol. 32 (7), pp. 913-923, 2014
 
Author: Chunming Li, all rights reserved
E-mail: li_chunming@hotmail.com
URL:  http://imagecomputing.org/~cmli/

This package includes two folders: MICO_2D and MICO_3D.

The folder MICO_2D includes the Matlab code for 2D MICO and some 2D test images. 

The folder MICO_3D includes the following files:
1. Matlab code for 3D MICO.
2. Two 3D images in Nifti format. 
3. Some Matlab code for reading and saving images in Nifti format.


Usage for 3D MICO
Input: 
      images in nifti format.
Output: 
      The results of segmentation and bias field correction in two files in nifti format: the       segmented image and the bias field corrected image are saved as nifti files with postfixes       "_seg" and "_bc", respectively. 

Note: For visualization, the code only shows the result for one of the 2D slices, although it performs 3D segmentation and bias field estimation/correction.
   





