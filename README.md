# MP3
Medical software for Processing multi-Parametric images Pipelines

**MP³** is an open source software aimed to support preclinical and clinical researchers for image processing studies. This Matlab toolbox offers a graphical interface that helps to process medical images. Here are the most important tools included in MP³ :
- Import medical images (e.g. DICOM, NIFTI, Bruker, …) and convert them to nifti/json format
- Visualize images (2D to 5D data)
- Process (e.g. registration, model fit, filters, …)
- Analyze (e.g. ROIs, scatterplots, …) medical images. 

Standard processes from SPM or FSL can be called as well as custom fonctions. Complex postprocessing pipelines can be created and applied to mutiparametric (e.g. MRI T1w + MRI T2w,…), multidimensional (e.g. spatial +temporal+echoes,…) or multimodal (e.g. MRI +CT, …) data. Once defined in a single patient, these processes can be stored and applied to a larger database using parallel architectures. MP³ was initially created in the Grenoble Institute of Neurosciences (France) to process multiparametric MRI protocols (anatomy+relaxometry+perfusion) in cohorts of +500 animals.

SUMMARY Poster ([ENGLISH](https://github.com/nifm-gin/MP3/blob/master/tools/Pictures/Poster_SFRMBM_Brossard_MP3_English.pdf)) ([FRENCH](https://github.com/nifm-gin/MP3/blob/master/tools/Pictures/Poster_SFRMBM_Brossard_MP3.pdf)): poster presented at the 4th Congress of the SFRMBM (Strasbourg).

You can find the article linked to MP3 here: https://www.frontiersin.org/articles/10.3389/fninf.2020.594799/full

Do not forget to cite our article if you used MP3 :)

   Brossard, C., Montigon, O., Boux, F., Delphin, A., Christen, T., Barbier, E. L., & Lemasson, B.
   MP3: Medical software for Processing multi-Parametric images Pipelines
   Frontiers in Neuroinformatics, 2020, 14, 53.
 
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


### User Guide: How to...
Please have a look at our tutorials available on our youtube channel : https://www.youtube.com/playlist?list=PL-Tj6Wc9aE9x7i6s-RLetvNE0isnEsFm7

*  [Import my data?](https://github.com/nifm-gin/MP3/wiki/User_guide_import_data)

*  [Manage my database?](https://github.com/nifm-gin/MP3/wiki/User_guide_manage_database)

*  [Visualize my data?](https://github.com/nifm-gin/MP3/wiki/User_guide_visualize_data)

*  [Create a pipeline and execute it?](https://github.com/nifm-gin/MP3/wiki/User_guide_create_execute_pipeline)

*  [Create a new module?](https://github.com/nifm-gin/MP3/wiki/User_guide_create_module)

*  [Download modules?](https://github.com/nifm-gin/MP3/wiki/User_guide_download_modules)


### Data test
The documentation found in the wiki page is based on this dataset : https://github.com/nifm-gin/Data_Test_MP3 (MP3 project)

The documentation on video is based on a dataset downloadable her : https://drive.google.com/drive/folders/1xWsNJf3e4auL_WWknvvGbk0UqjcfAY7A?usp=sharing (as Raw MRI data) or here https://drive.google.com/drive/folders/1MImE1AkOOKhOIXreA_jrIcJJloEXFJu1?usp=sharing (as a MP3 Project created with MP3 v1.3)

### Modules repository
We provide some modules available at https://github.com/nifm-gin/MP3_Modules_Repository

# Other
[Brouillon Module Coreg](https://github.com/nifm-gin/MP3/wiki/Brouillon_module_coreg)

