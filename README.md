# muVES
###
*Algorithm to analyze 2D and 3D images of microvascular networks (Politecnico di Milano)*

**Active contributors**: Alberto Rota, Martina Senesi, Adelaide Stucchi, Irene Venturelli
**Supervisor**: Luca Possenti
**Past contributors**: 
**Copyright**: Alberto Rota, Luca Possenti
**Mailto**: <alberto2.rota@mail.polimi.it> , <luca.possenti@polimi.it>

-------------------------------------------------------
Software requirements
-------------------------------------------------------
µVES has been tested on Matlab2020a.
Matlab2019a or further version required.

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


