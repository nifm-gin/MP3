# MP3
Medical software for Processing multi-Parametric images Pipelines

**MP³** is an open source software aimed to support preclinical and clinical researchers for image processing studies. This Matlab toolbox offers a graphical interface that helps to import (e.g. DICOM, NIFTII, Bruker, …), visualize (2D to 5D data), process (e.g. registration, model fit, filters, …) and analyze (e.g. ROIs, scatterplots, …) medical images. Standard processes from SPM or FSL can be called as well as custom fonctions. Complex postprocessing pipelines can be created and applied to mutiparametric (e.g. MRI T1w + MRI T2w,…), multidimensional (e.g. spatial +temporal+echoes,…) or multimodal (e.g. MRI +CT, …) data. Once defined in a single patient, these processes can be stored and applied to a larger database using parallel architectures. MP³ was initially created in the Grenoble Institute of Neurosciences (France) to process multiparametric MRI protocols (anatomy+relaxometry+perfusion) in cohorts of +500 animals.

SUMMARY Poster ([ENGLISH](https://github.com/nifm-gin/MP3/blob/master/tools/Pictures/Poster_SFRMBM_Brossard_MP3_English.pdf)) ([FRENCH](https://github.com/nifm-gin/MP3/blob/master/tools/Pictures/Poster_SFRMBM_Brossard_MP3.pdf)): poster presented at the 4th Congress of the SFRMBM (Strasbourg).

### Requirements
In order to fully enjoy MP³, you have to meet the following requirements:
* Matlab 2017b or higher
* Toolboxes (Mandatory : Image Processing Toolbox -- Recommended : Statistics and Machine Learning Toolbox ; Parallel Computing Toolbox)
* Java 8 or higher
* Data to process in Bruker / DICOM / Philips / Nifti format


### Download
You can download MP³ from our GIT repository. Then just add the downloaded file to your Matlab path and type `MP3` in Matlab command window. The graphical interface of the Viewer will then be displayed.


### Developers guide
In order to push your modifications to the community, please contact benjamin.lemasson@univ-grenoble-alpes.fr to get the developer status. Please create on your own branch from the dev one and then ask to merge it in the dev branch.


### Quickstart guide
*  [If you cannot wait !](https://github.com/nifm-gin/MP3/wiki/Quickstart_guide)


### User Guide: How to...

*  [Import my data?](https://github.com/nifm-gin/MP3/wiki/User_guide_import_data)

*  [Manage my database?](https://github.com/nifm-gin/MP3/wiki/User_guide_manage_database)

*  [Visualize my data?](https://github.com/nifm-gin/MP3/wiki/User_guide_visualize_data)

*  [Create a pipeline and execute it?](https://github.com/nifm-gin/MP3/wiki/User_guide_create_execute_pipeline)

*  [Create a new module?](https://github.com/nifm-gin/MP3/wiki/User_guide_create_module)



# Other
[Brouillon Module Coreg](https://github.com/nifm-gin/MP3/wiki/Brouillon_module_coreg)

[Fonctions à développer](https://github.com/nifm-gin/MP3/wiki/futur_features)

[Dev: nouveau module](https://github.com/nifm-gin/MP3/wiki/dev_new_module)
