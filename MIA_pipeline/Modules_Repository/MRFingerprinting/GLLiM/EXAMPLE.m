%%% This example compare the GLLiM algorithm using Lw=0 or Lw=1 on the
%%% Stanford face dataset (isomap.stanford.edu/datasets.html)
%%%
%%% The function gllim_face_test learn a mapping from face poses to images
%%% using all except 8 random test images. Two algorithms are tested,
%%% namely hGLLiM-0 and hGLLiM-1, with isotropic equal constraints on
%%% covariance matrices Sigma.
%%%
%%% Then, the pose of the test images are estimated using inverse mapping,
%%% and the images are reconstructed form estimated low-dimensional
%%% variables using forward mapping.
%%%
%%% -The first 2 figures show the log-likehood of the algorithms along
%%% iterations.
%%% -Figure 3 has 3 columns:
%%%   -the first one correspond to test images
%%%   -the second one to reconstructions using Lw=0
%%%   -the third one to reconstructions using Lw=1
%%% -Figure 4 shows the image reconstructed using different values
%%% of W, and Lw=1. Each row correspond to a test image, while the W
%%% value is increased from the leftmost to the rightmost column.

load data/face_data.mat;
gllim_face_test(images,poses);