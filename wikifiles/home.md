Welcome to the $\mu$VES Wiki! 

This page will guide you through:
* [Installation](#Installation)
* [Usage](#Usage)
* [Troubleshooting](#Troubleshooting)

A page with an [Example](https://github.com/alberto-rota/muVES/wiki/Example) is also available.

# Installation
Weather this repo is downloaded as a zipped package and extracted or it is cloned from Git, make sure that it gets located inside of a folder that inside of the MATLAB *PATH*. 
The root MATLAB folder `C:\Users\<user>\Documents\MATLAB` is already added to the *PATH*.

Add additional folders to the path is done in: `HOME > Environment > Set Path > Add with Subfolders`

### Version
The code has been developed and tested in the **MATLAB_2020a** release. You can report version incompatibilities by opening an [Issue](https://github.com/alberto-rota/muVES/issues) 

### Requirements
Other than the basic MATLAB package, the following Add-Ons are required:
* Image Processing Toolbox
* Curve Fitting Toolbox
* Computer Vision Toolbox
* Deep Learning Toolbox

Installation is possible via: `HOME > Environment > Set Path > Add with Subfolders`

# Usage
$\mu$VES is an automated tool for the topological and morphological analysis of 3D images of microvascular networks acquired from confocal microscopy. The supported file extensions are `.oib` and `.nd2`. A surrogate 2D version of $\mu$VES is also available and it supports `.bmp` image files.

### 1. Set up the analysis parameters 
$\mu$VES is highly flexible in terms of prioritizing accuracy or computational time. To best exploit such feature, the `muVES_settings.txt` configuration file can be personalized with the required parameters.

*Example:* By increasing the `Radius Precision` parameters, the radius of a vessel will be calculated by accounting more points along the longitudinal coordinate, granting a higher accuracy but requiring more computation.

A verbose description of each parameter is provided [in the file itself](https://github.com/alberto-rota/muVES/blob/main/muVES%20settings.txt).


### 2. Run the analysis
Once the analysis parameters have been set up, $\mu$VES is ready to run: launch it by opening the file `muVES_3D.m` and then click on the **Run** button in the editor tab (or press F5 as the default keyboard shortcut).

You will be prompted to a File Explorer window that you'll use to navigate through your files: select one of the files in the supported formats and click *Open*.

$\mu$VES is now running, and its progress can be monitored in the command window. *NOTE* that $\mu$VES **does NOT run as a background MATLAB process**, therefore the command window cannot be used for other tasks.

The computational time that $\mu$VES requires is dependent on:
* The specs of your machine (CPU and RAM speed above all)
* The parameters specified in the `muVES_settings.txt` configuration file
* The quality of the raw input image
* The topological and morphological complexity of the input image

### 3. Gather the output data
Once you see message
```
> Results saved at: <name-of-the-inpout-file>.mat

ans = 

struct with fields:
    info: ...
     raw: ...
    flat: ...
     ...: ...
```
in the *Command Window*, $\mu$VES will have finished its analysis. 

**Output files are saved in the same directory of your input files and with its same name, but different extensions**. For example, if you performed an analysis on

`C:\Users\<user>\myfiles\filename.oib`

you will find output files in the folder 

`C:\Users\<user>\myfiles`

namely:
* `filename.mat`: Contains a MATLAB variable named `mvn` with numerical results of the analysis, as the average length, radius, eccentricity, *etc*. It contains also the segmented and the skeletonized images of the network, as well as a list of its branches and their connectivity. The file [muVES_legend.txt](https://github.com/alberto-rota/muVES/blob/main/muVES%20legend.txt) is a detailed description of all the fields of the struct `mvn`
* `filename.bmp`: Contains a 2D projection of the 3D image
* `filename_radius.pts` and `filename_bifurcation.pts`: Input files for a 3D Poiseuille-Darcy Finite Elements simulation

### 4. Visualizing Results
The repository also includes the script `muVES_visualize.m`. This script produces appealing graphs and visual representations for the analysis performed by $\mu$VES.

In the same way as any other MATLAB script, it must be run with a click on the **Run** button in the editor tab (or with F5 as the default keyboard shortcut).










