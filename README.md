# µVES
###
*Algorithm to analyze the morphology and topology of 2D and 3D images of microvascular networks (Politecnico di Milano)*

**Active contributors**: Alberto Rota

**Supervisor**: Luca Possenti

**Copyright**: Alberto Rota, Luca Possenti

**Mailto**: <alberto1.rota@polimi.it> , <luca.possenti@polimi.it>

***
μVES is an algorithm that autonomously analyses 2D and 3D images of microvascular networks. It provides a complete morphological analysis with information on the vessel's length, radius, tortuosity and eccentricity, as well as on the network topology with connectivity, bifurcations and dead ends.

From the raw microscopy data, muVES outputs a the segmented network, its skeleton, it connectivity map and a table with all morphological and topological information regarding each vessel.

<table>
<tr>
<td align=center> 2D Projection </td>
<td align=center> Segmentation </td>
<td align=center> Skeletonization </td>
<td align=center> Connectivity </td>
</tr>
<tr>
<td> <img src=readme\projection.png width=300> </td>
<td> <img src=readme\segmentation.png width=300> </td>
<td> <img src=readme\skeleton.png width=300> </td>
<td> <img src=readme\graph.png width=300> </td>
<!-- <td> <img src=readme\orange.png width=300> </td> -->
</tr>
<!-- <td align=center> Radius map </td> -->
<td align=center> Length map </td>
<td align=center> Tortuosity map </td>
<td align=center> Eccentricity map </td>
<tr>
<td> <img src=readme\radius.png width=300 > </td>
<td> <img src=readme\lenght.png width=300> </td>
<td> <img src=readme\tortuosity.png width=300> </td>
<td> <img src=readme\eccentricity.png width=300> </td>
<!-- <td> <img src=readme\hists.png width=300> </td> -->
</tr>
</table>

***
## μVES 2D
μVES was designed to analyze 3D images of microvascular networks. However, its features have been adapted to work with 2D images as well. 

The 2D version of this algorithm outputs the same information as the 3D version, with the exception of:
- Z coordinates are never present
- Vessel Eccentricities cannot be calculated
- The lateral area is only estimated, assuming perfectly circular cross sections

***

## How to use μVES
A detailed illustrated guide on how to use μVES, including software requirements and dependencies, is provided at [this repo's wiki](https://github.com/alberto-rota/muVES/wiki)

<p align="left"> 
<a href="https://github.com/alberto-rota/muVES/wiki">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme/wikibutton_light.png">
  <img style="vertical-align:middle" alt="NEARLab" src="readme/wikibutton_dark.png" width="200" > 
</picture>
</a> </p>

***
## How to cite μVES
If you use μVES in your research, please cite the following paper:

```
@article{rota2023muves,
author = {Rota, Alberto and Possenti, Luca and Offeddu, Giovanni S. and Senesi, Martina and Stucchi, Adelaide and Venturelli, Irene and Rancati, Tiziana and Zunino, Paolo and Kamm, Roger D. and Costantino, Maria Laura},
title = {A three-dimensional method for morphological analysis and flow velocity estimation in microvasculature on-a-chip},
journal = {Bioengineering \& Translational Medicine},
volume = {n/a},
number = {n/a},
pages = {e10557},
keywords = {3D computational analysis, deep learning, network morphology, segmentation, vasculature-on-a-chip},
doi = {https://doi.org/10.1002/btm2.10557},
url = {https://aiche.onlinelibrary.wiley.com/doi/abs/10.1002/btm2.10557},
}

```

<!-- -------------------------------------------------------
Software requirements
-------------------------------------------------------
µVES has been tested on Matlab2020a.
Matlab2020a or further version required, with *Image Processing Toolbox*, *Curve Fitting
Toolbox*, *Computer Vision Toolbox* and *Image Processing Toolbox*

-------------------------------------------------------
Setting parameters
-------------------------------------------------------
Prior to running the script, please specify all the analysis settings in the file 
"muVES settings.txt", most importantly check the pixel resolution in micrometers/pixel.

-------------------------------------------------------
Run the analysis
-------------------------------------------------------
For a **single image** analysis open and run the script "muVES_3D". A file selection window
will open, from which you will select a file of supported type [.oib, .n2d]. When the analysis
is complete, the results will be saved in the same folder of the selected files: these are
a .mat file with all the results (its contents are covered in the file "muVES legend.txt"),
a flat 2D projection, and the .pts files to be used for the CFD analysis.

"muVES_3D_process_batch" is also available: after running it, a **folder selection window** will 
open, from which you will select a folder of interest. All files of supported type contained in
the folder will be analysed, and a "batch.xlsx" file will be provided, with the "branchdata" 
table from each file analyzed in its own sheet. 

N.B.: Both of these files are available for the **2D analysis**. The supported filetypes are [.oib, .nd2,
.bmp], but 3D data [.oib, .nd2] will be projected prior to the analysis

-------------------------------------------------------
Visualizing results
-------------------------------------------------------
The script "muVES_visualize" is a visualization tool for all the graphic informations obtained 
from the analysis. Specify (in the first lines of the script) which plots/graphs/renderings you
want to be shown and then run the script. From the file selection window that opens, please 
select the .mat file that you want to visualize.

-------------------------------------------------------
Analyzing flowrate
-------------------------------------------------------
Flowrate in each vessels can be computed by defining pressure BC (inlet pressure). 
Outlet pressure is equal to 0. Simmetry condition is defined as null flow rate.
Define the BC in the script, and run it. Choose the .mat file resulting from the previous anlysis. -->
