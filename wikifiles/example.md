# $\mu$VES - An Example

This an illustrated step-by-step guide in performing an analysis with $\mu$VES.

### 1. Setting the analysis parameters
1. Open `muVES_settings.txt`
2. Change the values in lines *1-11*. For example, if your image has a high resolution and wouldn't suffer from a downsampling operation, set `Downsampling Factor` to 2.

![settingsimg](https://github.com/alberto-rota/muVES/wikifiles/settings.png)

### 2. Running the algorithm        
1. Open `muVES_3D.m`
2. Click the `Run` button
3. Once the file explorer window is open, sselect the file you want to analyze. `test.oib` is provided in this repository in the *Test Images* folder
4. Click `Open`
5. Monitor the MATLAB Command Window to see the progress

![runningimg](https://github.com/alberto-rota/muVES/wikifiles/running.png)

6. When $\mu$VES has finished running, alongside your input file (in the same directory) you'll see the output files generated.
7. The Command Windows will also display the content of `test.mat`
   
![outputfilesimg](https://github.com/alberto-rota/muVES/wikifiles/outputfiles.png)

### 3. Visualizing results