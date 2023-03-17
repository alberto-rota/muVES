# MVN Legend
* info: 		General information about the network
   * chrono: 		Computational time used for all the steps of the analysis
   * hashtable:		Java Hashtable (as a string) with details from microscopy
   * name:           Path of the input data
   * pxdensity: 		Space resolution in the 3 dimensions (in micrometers/pixel)
   * downfactor: 	Downsampling factor
   * pts_per_branch:	Number of voxel for each branch used in the CFD simulation
   * ragPrecision: 	Number of points for each branch used to calculate the radius
   * lungThr:		Branches shorter that this value (in voxel) are removed
   * ptsOK:          True if the .pts is written correctly
* raw: 		Input 3D volume 
* flat: 		2D projection of all the layers
* bw: 		Binarized volume obtained after segmentation
* skel: 		Topological results
   * sk: 		Binary matrix with the skeleton
   * bp:			Binary matrix with the branchpoints
   * ep:			Binary matrix with the endpoints
   * graph:		Graph of the network
* branchdata: 	Table with the values of all the branches
   * Num:		Branch number (it is also its identifier)
   * xPath:		X coordinates of the voxels in the branch
   * yPath:		Y coordinates of the voxels in the branch
   * zPath:		Z coordinates of the voxels in the branch
   * From:		Starting point of the branch
   * CatFr: 		Category of the starting point [INT/MIX/DIR] 
   * To: 		Ending point of the branch
   * CatTo: 		Category of the ending point [INT/MIX/DIR] 
   * Interp: 		Interpolation parameters for the branch (piecewise polynomial)
   * xInt:		X coordinates of the voxels used for the CFD
   * yInt:		Y coordinates of the voxels used for the CFD
   * zInt:		Z coordinates of the voxels used for the CFD
   * subN:		Subnetwork to which the branch belongsù
   * Len:		Length of the branch
   * Tort:		Tortuosity of the branch
   * Rad:		Radius of the branch
   * Alat:		Lateral surface area (as a cylinder) of the branch
   * Eccent:		Eccentricity of the branch
   * Res: 		Approximation of the hydraulic resistance
   * isGood:		Flag for the goodness of the branch. TRUE if length > 3*radius 
* mRad: 		Mean value of the radius of the branches
* mLen: 		Mean value of the lenght of the branches
* mTort: 	Mean value of the tortuosity of the branches
* mEcc: 		Mean value of the eccentricity of the branches
* volFrac:	% of white voxels
* approxAlat:	Lateral surface area with the cilindical approximation
* realAlat:	Lateral surface area with counting the number of voxels
* S_over_V: 	Ratio SurfaceArea/Volume
* numSubN: 	Number of disconnected subnetworks 
* pcr:		PCR