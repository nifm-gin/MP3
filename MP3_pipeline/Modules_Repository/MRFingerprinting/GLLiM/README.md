
Files from https://github.com/Chutlhu/GLLiM initial contribution

# GLLiM
a flexible Matlab toolbox for Gaussian Locally Linear Mapping 

The GLLiM toolbox provides a set of Matlab functions allowing to learn a relationship between two spaces. It implements the  hGLLiM algorithm described in detail in: A. Deleforge, F. Forbes, and R. Horaud. High-Dimensional Regression with Gaussian Mixtures and Partially-Latent Response Variables. Statistics and Computing. 2014. The article is available online on arXiv at http://arxiv.org/abs/1308.2302.

The toolbox also contains an example of application (EXAMPLE.m) using the freely available Stanford dataset for faces: http://isomap.stanford.edu/datasets.html.


# TODO:
- Verify that the normalization term in is corrected ( sqrt() )
- gllim_inverse_dens_modified.m and gllim_direct_map.m always useful ?
- cov matrix parallelizable ? (rename variables)
