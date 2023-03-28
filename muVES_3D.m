% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto1.rota@polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@mail.polimi.it
% ========================================================================%
% SPECIFY THE SETTINGS IN THE "muVES settings.txt" FILE 
% ======================================================================= %
%% SETUP AND IMAGE LOADING
function mvn = muVES_3D(varargin)
% ======================================================================= %
% REQUIREMENTS: 
req = {'Image Processing Toolbox', 'Curve Fitting Toolbox', ...
    'Computer Vision Toolbox','Deep Learning Toolbox'};
% ======================================================================= %
if nargin == 0
% Specify the path to the file. The extension must be separately specified
% in 'extension'. If those fields are left empty ('') or if the path is not
% correct, a window for file selection will be opened.
% Example:      pathtoimg = 'C:\...\...\myfolder\myfile';
%               extension = '.oib';
pathtoimg = "";
extension = ".oib";
else
    splpath = strsplit(varargin{1},".");
    pathtoimg = splpath(1);
    extension = strcat(".",splpath(2));
end
fid = fopen("muVES settings.txt");
% ========================================================================%
% 'pxdens' contains the space resolution of the microscope in the 3
% dimensions, specified in micrometers/pixel. The algorithm interpolates
% the volume so that the resolution in the Z direction is the same to the
% resoution in the X and Y direction
xds = strsplit(fgetl(fid),":");
yds = strsplit(fgetl(fid),":");
zds = strsplit(fgetl(fid),":");
pxdens = [str2double(xds{2}) str2double(yds{2}) str2double(zds{2})];
% [POSITIVE REAL 3x1 VECTOR]
% ========================================================================%
% The number of voxels in the 3D matrix gets reduced of a factor
% 2^downfactor, keeping 1-every-downfactor voxels in each dimensions.
% 'downfactor = 1' keeps the matrix unaltered, 'downfactor = 2' keeps 1
% voxel every 2.
pds = strsplit(fgetl(fid),":");
downfactor = str2double(pds{2}); % [POSITIVE INTEGER SCALAR]
% [INCREASE for PERFORMANCE, DECREASE for ACCURACY]
% By doubling the downfactor, the number of operations is approximately 
% reduced by 8 times and with it, the computational time required.
% ========================================================================%
% 'smoothing_repeat' indicates how many consecutive times a smoothing
% operation is performed on the binary matrix.
% The higher this value is, the more times the matrix is ​​"refined", making it
% clearer and decreasing the skeletonization errors; at the same time some
% information is lost(especially on thin vessels)
pds = strsplit(fgetl(fid),":");
smoothing_repeat = str2double(pds{2});  % [POSITIVE INTEGER SCALAR]
% ========================================================================%
% 'rad_precision' is the number of equally spaced sections obtained on
% each branch for the calculation of the radius and lateral area
pds = strsplit(fgetl(fid),":");
rad_precision = str2double(pds{2});   % [POSITIVE INTEGER SCALAR]
% [INCREASE for ACCURACY, DECREASE for PERFORMACE]
% ========================================================================%
% Set 'print_pts' to true or false if the .pts file has to be
% produced
pds = strsplit(fgetl(fid),":");
print_pts = str2double(pds{2});  % [LOGICAL SCALAR]
% 'pts_per_branch' (INTEGER) is the number of voxels in eanch branch in
% the skeleton provided to the CFD simulation
pds = strsplit(fgetl(fid),":");
pts_per_branch = str2double(pds{2});  % [POSITIVE INTEGER SCALAR]
% 'edge_border' [FLOAT from 0 to 1] is ratio between the edge width and the
% total width of the image. 'edge_border = 0.05' means that 5% of the total
% width is considered edge (used to apply the boundary conditions)
pds = strsplit(fgetl(fid),":");
edge_border = str2double(pds{2}); % [REAL SCALAR (in the range 0-1)]
% ========================================================================%
% 'lung_thr' is the minimum length for a branch to be
% kept in the skeleton
pds = strsplit(fgetl(fid),":");
lung_thr = str2double(pds{2}); % [POSITIVE REAL SCALAR]
% ========================================================================%
% All the results are stored in the struct 'mvn', accessible in 
% the 'filename.mat' created in the same folder as the input.
% --> See 'MVN legend.txt' for all the details of the values in 'mvn'. 
% To specify something to append to the .mat with the results, add it to
% 'append_to_path', otherwise leave it as "".
pds = strsplit(fgetl(fid),":");
append_to_path = num2str(pds{2}); 
% ========================================================================%
fclose(fid);
if ~check_req(req)
    disp('Requirements:');
    disp(string(req'));
    error('ERROR: You are missing some toolboxes that are required for the script');
end
try
%---------------------------------------------------- ---------------------%
% Melissa Linkert, Curtis T. Rueden, Chris Allan, Jean-Marie Burel, 
% Will Moore, Andrew Patterson, Brian Loranger, Josh Moore, Carlos Neves,
% Donald MacDonald, Aleksandra Tarkowska, Caitlin Sticco, Emma Hill, Mike
% Rossner, Kevin W. Eliceiri, and Jason R. Swedlow (2010) 
% "Metadata matters: access to image data in the real world." 
% The Journal of Cell Biology 189(5), 777-782. doi: 10.1083/jcb.201004104
%-------------------------------------------------------------------------%
    img_3d = bfOpen3DVolume(char(strcat(pathtoimg,extension)));
    vol = img_3d{1,1}{1,1};
    splpath = strsplit(img_3d{1}{2},".");
    pathtoimg = splpath(1);
    extension = strcat(".",splpath(2));
catch
    error('Invalid path or unreadable file');   
end
% ========================================================================%
% SPECIFY HERE THE ADDITIONAL OPERATION THAT HAVE TO BE PERFORMED FOR
% SPECIFIC FILE EXTENSION. The final data must be named 'vol'
switch extension
    case {'.oib'}
        
    case {'.nd2'}
        vol = vol(:,:,1:round(end/4));
end
tic
pathtoimg = char(pathtoimg);
disp(strcat("> Processing: ",string(pathtoimg),extension));
% ========================================================================%
% % WIP - CONTRAST ENHANCING
% vol = mat2gray(vol);
% % Contrast enhancing. We keep black voxels dark while increasing the
% % brightness of gray and white pixels
% vol = imadjustn(vol,[],[],0.1);
% vol = vol.*1000;
% ========================================================================%
% Chrono initialization
mvn.info.chrono = table("File reading and initialization",toc,'VariableNames',...
    {'Phase','Duration'});
mvn.info.hashtable = string(img_3d{1,2});
mvn.raw = vol;
clear img_3d splpath path  xds yds zds pds fid

%% 2D PROJECTION
tic
disp("> Creating 2D projection ");
% 2D projection
flattened = max(vol, [],3);
mvn.flat = flattened;
% Saving the 2D projection as 'filename_flat.bmp'
truec = cat(3, flattened,flattened,flattened);
imwrite(mat2gray(truec),[pathtoimg '.bmp']);

try
    mvn.info.chrono = [mvn.info.chrono; {"2D Flattening ",toc}];
catch
end
clear truec ca flattened

%% DOWNSAMPLING AND REPROPORTIONING
tic
disp("> Downsampling and adjusting the voxel dimensions ");
% Downsampling operation, selecting one-every-downfactor voxels along the
% tree dimensions of the 3D matrix
V = vol(1:downfactor:end, 1:downfactor:end, 1:downfactor:end);

% Dimension of the downsampled volume, Height, Witdth and Layers
h = size(V, 1);
w = size(V, 2);
l = size(V, 3);

% Creating the meshgrid
[x,y,z] = meshgrid(1:w,1:h,1:l);

% Voxel dimension adjustment. If the pixel density is different in the
% vertical direction, the voxel represent a parallelepiped, and are
% converted to a cube with a linear interpolation
if pxdens(3)~=pxdens(1)
    newvol = [];
    for slh = linspace(1,l,(l-1)*pxdens(3)/pxdens(1))
        s = slice(x,y,z,double(V),[],[],slh,'linear');
        newvol = cat(3,newvol,(s.CData)');
    end
    s.Visible = false;
    vol = smooth3(newvol,'gaussian');
end

try
    mvn.info.chrono = [mvn.info.chrono; ...
        {"Downsampling and voxel dimension adjustment",toc}];
catch
end

close
clear V s slh newvol

%% SEGMENTATION
% If the "Current Folder" contains "mVN_DLN.mat" the segmantation is
% performed with a faster Deep Learning technique. Otherwise, ActiveContour
% is used. To pureposely use ActiveContour, instead, please comment out
% "load mVN_DLN" in the following lines of code
tic

% INITIALIZING SEGMENTATION WITH DEEP LEARNING
load mVN_DLN;
if exist('mVN_DLN','var')   
    disp("> Performing segmentation - Deep Learning ");
    
    % The image is scaled to 'rescfact' to be best interfaced with the
    % neural network. It will then be re-scaledo its original size after
    % segmentation (must be a multiple of 128x128x8)
    % CHANGE THIS RESCALING ACCORDING TO COMPUTATIONAL POWER
    rescfact = [384 384 16];
    vres = imresize3(vol,rescfact);
    
    % Segmenting quadrants of dimension 128x128x8. The final segmentation
    % is the re-alignment of all these quadrants
    resultres = zeros(rescfact);
    for i=0:128:(rescfact(1)-128)
        for j=0:128:(rescfact(2)-128)
            for k=0:8:(rescfact(3)-8)
                resultres((i+1):(i+128), (j+1):(j+128), (k+1):(k+8)) = ...
                    semanticseg(vres((i+1):(i+128), (j+1):(j+128), (k+1):(k+8)),mVN_DLN);
            end
        end
    end
    % Rescaling the image up to its original size
    resultres = imresize3(resultres,size(vol));
    bw = resultres  >= 1.5;
    
    % Smoothing edges
    for i=1:smoothing_repeat
        bw = smooth3(bw);
    end
    bw = bw > 0.5;
    
    % Checking a possible missegmentation with DL. If such, segmentation is
    % performed again with AC
    if nnz(bw)/numel(bw) < 0.05
        clear bw;
        warning("Bad segmentation with DL. Re-Segmenting with AC");
        disp("> Reperforming segmentation - ActiveContour ");
        mask = zeros(size(vol));
        
        for i=1:smoothing_repeat
            vol = smooth3(vol);
        end
        
        for seedLevel = 1:size(vol,3)
            % The mask is initialized by comparing all pixels in the brighter layer
            % with a threshold value, obtaining a binary image that is the
            % 2D segmentation (inaccurate) of this layer.
            seed = vol(:,:,seedLevel) > mean(vol,'all');
            % The mask gets thickened by 10 pixels and then furtherly refined
            % with a moving average 2D filter of arbitrary size of 10x10.
            seed_thick = bwmorph(seed, 'thicken', 10);
            N = 3;
            kernel = ones(N, N) / N^2;
            seed_thick = conv2(double(seed_thick), kernel, 'same');
            
            % The result of the moving average is of type double.In order to
            % reconvert it to type 'logical' it must be compared with a
            % threshold value that depends on convolution kernel size
            seed_thick = seed_thick > N/10;
            mask(:,:,seedLevel) = seed_thick;
        end
        
        try
            mvn.info.chrono = [mvn.info.chrono; {"Mask Creation",toc}];
        catch
        end
        
        % Actual activecontour segmentation
%-------------------------------------------------------------------------%
% Chan T., Vese L. (1999)
% "An Active Contour Model without Edges"
% In: Nielsen M., Johansen P., Olsen O.F., Weickert J. (eds) Scale-Space
% % Theories in Computer Vision. Scale-Space 1999. Lecture Notes in Computer
% Science, vol 1682. Springer, Berlin, Heidelberg.
% https://doi.org/10.1007/3-540-48236-9_13
%-------------------------------------------------------------------------%
        bw = activecontour(vol,mask,100);
        
        try
            mvn.info.chrono = [mvn.info.chrono; {"ReSegmentation - ActiveContour",toc}];
        catch
        end
    else
        try
            mvn.info.chrono = [mvn.info.chrono; {"Segmentation - Deep Learning",toc}];
        catch
        end
    end
    
else
    % SEGMENTATION WITH ACTIVECONTOUR
    % The mask used is the binary image of the layer with the overall higher
    % brightness level. The layer considered is identified by
    % 'seedLevel'. This operation is done only on the smoothed volume
    
    tic
    disp("-DL architecture not found");
    disp("> Performing segmentation - ActiveContour ");
    mask = zeros(size(vol));
    
    for i=1:smoothing_repeat
        vol = smooth3(vol);
    end
    
    for seedLevel = 1:size(vol,3)
        % The mask is initialized by comparing all pixels in the brighter layer
        % with a threshold value, obtaining a binary image that is the
        % 2D segmentation (inaccurate) of this layer.
        seed = vol(:,:,seedLevel) > mean(vol,'all');
        % The mask gets thickened by 10 pixels and then furtherly refined
        % with a moving average 2D filter of arbitrary size of 10x10.
        seed_thick = bwmorph(seed, 'thicken', 10);
        N = 3;
        kernel = ones(N, N) / N^2;
        seed_thick = conv2(double(seed_thick), kernel, 'same');
        
        % The result of the moving average is of type double.In order to
        % reconvert it to type 'logical' it must be compared with a
        % threshold value that depends on convolution kernel size
        seed_thick = seed_thick > N/10;
        mask(:,:,seedLevel) = seed_thick;
    end
    
    try
        mvn.info.chrono = [mvn.info.chrono; {"Mask Creation",toc}];
    catch
    end
    
    % Actual activecontour segmentation
%-------------------------------------------------------------------------%
% Chan T., Vese L. (1999) 
% "An Active Contour Model without Edges"
% In: Nielsen M., Johansen P., Olsen O.F., Weickert J. (eds) Scale-Space 
% % Theories in Computer Vision. Scale-Space 1999. Lecture Notes in Computer 
% Science, vol 1682. Springer, Berlin, Heidelberg. 
% https://doi.org/10.1007/3-540-48236-9_13
%-------------------------------------------------------------------------%
    bw = activecontour(vol,mask,100);
    
    try
        mvn.info.chrono = [mvn.info.chrono; {"Segmentation - ActiveContour",toc}];
    catch
    end
end

clear ca vhandle mask kernel mVN_DLN N rescfact resultres vres seedLevel ...
    seed_thick seed i j k

%% ALIGNMENT ON THE HORIZONTAL PLANE
tic
disp("> Aligning on the horizonal plane ");
[h,w,l] = size(bw);
% Computing the plane approximating the whole network
[x_bw,y_bw,z_bw]=ind2sub([h,w,l],find(bw));
f = fit([x_bw y_bw],z_bw,'poly11');

% Coefficients of the interpolating plane, in the form z = q + mx*x + my*y
q = f.p00; mx = f.p10; my = f.p01;

% The alignment value at every point is the difference between the height
% of the point and the maximum value of 'alignmat'. Each column vector 
% is shifted by a quantity equal to that of 'alignmat'
% alignmat = alignmat-max(max(alignmat));

bwalign = zeros(w,h,l);
for i=1:h
    for j=1:w
        if any(bw(i,j,:))
            newh = l+abs(round(-(q+mx*i+my*j)));
            bwalign(i,j,1:newh) = addshift(squeeze(bw(i,j,:)),round((q+mx*i+my*j)),1);
        end
    end
end

bw = logical(bwalign);
mvn.bw = bw;

try
    mvn.info.chrono = [mvn.info.chrono; {"Alignment on the horizontal plane ",toc}];
catch
end

clear newh bwalign alignmat i j f x_bw y_bw x_bw z_bw q mx my

%% SKELETONIZATION - SKELETON3D
tic
disp("> Computing the skeleton ");
% Computing the skeleton with 'Skeleton3D'
%-------------------------------------------------------------------------%
% Copyright (c) 2016, Philip Kollmannsberger
% All rights reserved.

% Ta-Chih Lee, Rangasami L. Kashyap and Chong-Nam Chu
% "Building skeleton models via 3-D medial surface/axis thinning 
% algorithms."
% Computer Vision, Graphics, and Image Processing, 56(6):462–478, 1994.

% Kerschnitzki, Kollmannsberger et al.,
% "Architecture of the osteocyte network correlates with bone material 
% quality."
% Journal of Bone and Mineral Research, 28(8):1837-1845, 2013.
%-------------------------------------------------------------------------%
sk = Skeleton3D(bw);
% For some reason Skeleton3D swaps the X and Y dimensions
sk = permute(sk, [2 1 3]);

mvn.skel.sk = sk;

try
    mvn.info.chrono = [mvn.info.chrono; {"Skeletonization",toc}];
catch
end

clear ca

%% COMPUTING BRANCHPOINTS AND ENDPOINTS
tic
disp("> Finding branchpoints ");
bpoints = bwmorph3(sk,'branchpoints');
epoints = bwmorph3(sk,'endpoints');
mvn.skel.bp = bpoints;
mvn.skel.ep = epoints;

try
    mvn.info.chrono = [mvn.info.chrono; {"Finding Endpoints and Branchpoints",toc}];
catch
end

clear bpoints epoints

%% INTERPOLATION AND NODE CLASSIFICATION
tic
disp("> Interpolating the branches ");
% Branches are split with the 'Skel2Graph3D' function. The information
% on the branches is contained in the vector of struct 'link'. The
% information about the various nodes is contained in the vector of struct
% 'node'.
%-------------------------------------------------------------------------%
% Copyright (c) 2016, Philip Kollmannsberger
% All rights reserved.

% Kerschnitzki, Kollmannsberger et al.,
% "Architecture of the osteocyte network correlates with bone material 
% quality."
% Journal of Bone and Mineral Research, 28(8):1837-1845, 2013.
%-------------------------------------------------------------------------%
[adj, node, link] = Skel2Graph3D(sk,lung_thr);

[h,w,l] = size(mvn.bw);

% Il numero di branches è dato dalla quantià di elementi di 'link'
tot_branches = numel(link);

% Viene creata una table 'vessel data' che contiene le seguenti
% informazioni:
%   -'Number' identifies the branch
%   -'xPath','yPath' and 'zPath' are the coordinates of all the points in
%     the branch
%   -'From' contain the coordinates of the starting point
%   -'To' contain the coordinates of the final point
%   -'CatFr' and 'CatTo' contain the categories ('INT', 'MIX', 'DIR')
%     of the branchpoint
%   -'Interp' contains the interpolation parameters
%   -'xInt','yInt' and 'zInt' are the coordinates of the points used in
%     the .pts file.

branchdata = table(0,{0},{0},{0},[0,0,0],{0},[0,0,0],{0},{0},...
    zeros(1,pts_per_branch-2),zeros(1,pts_per_branch-2),zeros(1,pts_per_branch-2),...
    'VariableNames', {'Num','xPath','yPath','zPath','From','CatFr','To','CatTo',...
    'Interp','xInt','yInt','zInt'});

% For each branch, the coordinates of the discrete points in the skeleton
% are extracted and then added to 'branchdata'
tic
for b=1:tot_branches
    path_idx = link(b).point;
    [path_x, path_y, path_z] = ind2sub([h,w,l],path_idx);
    interp = cscvn([path_x; path_y; path_z]);
    x_from = node(link(b).n1).comx;
    y_from = node(link(b).n1).comy;
    z_from = node(link(b).n1).comz;
    x_to = node(link(b).n2).comx;
    y_to = node(link(b).n2).comy;
    z_to = node(link(b).n2).comz;
    from = [x_from,y_from,z_from];
    to = [x_to,y_to,z_to];
    toAdd = {b, path_x, path_y, path_z, from,"", to, "",interp, ...
        zeros(1,pts_per_branch-2),zeros(1,pts_per_branch-2),...
        zeros(1,pts_per_branch-2)};
    branchdata = cat(1, branchdata, toAdd);
end
branchdata(1,:) = [];

% GRAPH CONVERSION
% From the adjacency matrix obtained from 'Skel2Graph3D', the graph of the
% network is extracted
tic
G = graph(adj);
for n = 1:numel(node)
    G.Nodes.x(n) = node(n).comx;
    G.Nodes.y(n) = node(n).comy;
    G.Nodes.z(n) = node(n).comz;
    G.Nodes.subN(n) = 0;
end

% The 'floodgraph' function distinguishes different disconnected sub-regions in the
% graph. Each node on the network is assigned the value 'subN', which
% identifies the sunetworks that it belongs to. See the description
% of the floodgraph function in the last section.

numsn = 1;
while any(G.Nodes.subN == 0)
    non_lab = find(G.Nodes.subN==0);
    G = floodgraph(G,non_lab(1),numsn);
    numsn=numsn+1;
end
mvn.skel.graph = G;

try
    mvn.info.chrono = [mvn.info.chrono; {"Interpolation + Graph",toc}];
catch
end

if print_pts
    tic
% For each branch, a vector of 'pts_per_brunch' points is created
% independently from its length. 
for i=1:tot_branches
    tmax = max(branchdata.Interp{i}.breaks);
    interv = linspace(tmax/(pts_per_branch-2),tmax-tmax/(pts_per_branch-2),...
        (pts_per_branch-2));
    for j = 1:pts_per_branch-2
        p = ppval(branchdata.Interp{i}, interv);
        branchdata.xInt(i,:) = p(1,:);
        branchdata.yInt(i,:) = p(2,:);
        branchdata.zInt(i,:) = p(3,:);
    end
end

% Etracting start and end point coordinate from each branch
s = branchdata.From;
e = branchdata.To;

% Each node of the skeleton is categorized with the appropriate condition
% required by the CFD input protocol:
% -"INT" for branchpoints
% -"DIR" for endpoints near edges, with a pressure value
% dependent on the nearest edge (1 to 4)
% -"MIX" for endpoints far form the edges

G.Nodes.deg = degree(G);

    % THIS FOR CYCLE ASSIGNS THE CONDITIONS TO THE BRANCHPOINTS [MIX-DIR-INT]
    for i = 1: tot_branches
        % Identifying branchpoints in the graph (order > 1)
        if G.Nodes.deg(all(branchdata.From(i,:) == table2array(G.Nodes(:,1:3)),2)) == 1
            flag_bpST = 0;
        elseif G.Nodes.deg(all(branchdata.From(i,:) == table2array(G.Nodes(:,1:3)),2)) > 1
            flag_bpST = 1;
        end
        
        if flag_bpST == 0
            if s(i,1)/h > (1-edge_border)
                branchdata.CatFr(i) = "DIR 1";
            elseif s(i,2)/w > (1-edge_border)
                branchdata.CatFr(i) = "DIR 2";
            elseif s(i,1)/h < edge_border
                branchdata.CatFr(i) = "DIR 3";
            elseif s(i,2)/w < edge_border
                branchdata.CatFr(i) = "DIR 4";
            else
                branchdata.CatFr(i) = "MIX";
            end
        else
            branchdata.CatFr(i) = "INT";
        end
        
        if G.Nodes.deg(all(branchdata.To(i,:) == table2array(G.Nodes(:,1:3)),2)) == 1
            flag_bpEN = 0;
        elseif G.Nodes.deg(all(branchdata.To(i,:) == table2array(G.Nodes(:,1:3)),2)) > 1
            flag_bpEN = 1;
        end
        
        if flag_bpEN == 0
            if e(i,1)/h > (1-edge_border)
                branchdata.CatTo(i) = "DIR 1";
            elseif e(i,2)/w > (1-edge_border)
                branchdata.CatTo(i) = "DIR 2";
            elseif e(i,1)/h < edge_border
                branchdata.CatTo(i) = "DIR 3";
            elseif e(i,2)/w < edge_border
                branchdata.CatTo(i) = "DIR 4";
            else
                branchdata.CatTo(i) = "MIX";
            end
        else
            branchdata.CatTo(i) = "INT";
        end
    end
    
    branchdata.CatFr = categorical(branchdata.CatFr);
    branchdata.CatTo = categorical(branchdata.CatTo);
    
    for i=1:tot_branches
        idx = all(branchdata.From(i,:) == [G.Nodes.x G.Nodes.y G.Nodes.z],2);
        branchdata.subN(i) = G.Nodes.subN(idx);
    end
    
    for i=1:tot_branches
        if branchdata.CatFr(i) == "DIR 1" || branchdata.CatTo(i) == "DIR 1" || ...
                branchdata.CatFr(i) == "DIR 2" || branchdata.CatTo(i) == "DIR 2" ||...
                branchdata.CatFr(i) == "DIR 3" || branchdata.CatTo(i) == "DIR 3" ||...
                branchdata.CatFr(i) == "DIR 4" || branchdata.CatTo(i) == "DIR 4"
            branchdata([1,i],:) = branchdata([i,1],:);
            break;
        end
    end
    
    if branchdata.CatTo(1) == "DIR 1" || branchdata.CatTo(1) == "DIR 2" || ...
            branchdata.CatTo(1) == "DIR 3" || branchdata.CatTo(1) == "DIR 4"
        appPoint = branchdata.From(1,:);
        branchdata.From(1,:) = branchdata.To(1,:);
        branchdata.To(1,:) = appPoint;
        appCat = branchdata.CatFr(1);
        branchdata.CatFr(1) = branchdata.CatTo(1);
        branchdata.CatTo(1) = appCat;
        branchdata.xInt(1,:) = flip(branchdata.xInt(1,:));
        branchdata.yInt(1,:) = flip(branchdata.yInt(1,:));
        branchdata.zInt(1,:) = flip(branchdata.zInt(1,:));
    end
    
    for i=2:tot_branches
        if branchdata.CatFr(i) == "DIR 1" || branchdata.CatFr(i) == "DIR 2" ||...
                branchdata.CatFr(i) == "DIR 3" || branchdata.CatFr(i) == "DIR 4"
            appPoint = branchdata.From(i,:);
            branchdata.From(i,:) = branchdata.To(i,:);
            branchdata.To(i,:) = appPoint;
            appCat = branchdata.CatFr(i);
            branchdata.CatFr(i) = branchdata.CatTo(i);
            branchdata.CatTo(i) = appCat;
            branchdata.xInt(1,:) = flip(branchdata.xInt(1,:));
            branchdata.yInt(1,:) = flip(branchdata.yInt(1,:));
            branchdata.zInt(1,:) = flip(branchdata.zInt(1,:));
        end
    end
    
    for i = 1:max(branchdata.subN)
        cond = branchdata.subN == i;
        idxx = find(cond);
        if all(branchdata.CatFr(cond) == "MIX" & branchdata.CatTo(cond) == "INT"| ...
                branchdata.CatFr(cond) == "INT" & branchdata.CatTo(cond) == "MIX"| ...
                branchdata.CatFr(cond) == "MIX" & branchdata.CatTo(cond) == "MIX"| ...
                branchdata.CatFr(cond) == "INT" & branchdata.CatTo(cond) == "INT") ...
                && numel(idxx) > 0
%             changeable = find(branchdata.CatFr == "MIX" & cond);
            branchdata(cond,:) = [];
        end
    end
    tot_branches = size(branchdata,1);
end
mvn.branchdata = branchdata;
try
    mvn.info.chrono = [mvn.info.chrono; {"Node classification",toc}];
catch
end

% Updating the discrete skeleton to the one maintained with the graph
% conversion. Short vessels have been eliminated 
sk = zeros(h,w,l);
for i=1:tot_branches
    for j=1:numel(branchdata.xPath{i})
        sk(branchdata.xPath{i}(j),branchdata.yPath{i}(j),branchdata.zPath{i}(j)) = 1;
    end
end


clear adj appCat appPoint b cond e flag_bpEN flag_bpST from G i idx idxx ...
    interp interv j link node p path_idx path_x path_y path_z s tmax ...
    to toAdd ncl non_lab x_from x_to y_from y_to z_from z_to

%% MORPHOLOGICAL MEASUREMENTS
tic
disp("> Morphological measurements ");
displine = 0;
errors = 0;
for b=1:tot_branches
    % LENGTH
    xyzpath = [branchdata.xPath{b};branchdata.yPath{b};branchdata.zPath{b}];
    branchdata.Len(b) = sum(sqrt(sum(diff(xyzpath,1,2).^2)));
    
    % TORTUOSITY [as T = Dist/L]
    d = sqrt(sum((branchdata.From(b,:)-branchdata.To(b,:)).^2));
    branchdata.Tort(b) = constrain(branchdata.Len(b)/d,[1 Inf]);
    
    % TORTUOSITY [new approach as variance of the mean angle]
    xyzdir = diff(xyzpath,1,2);
    theta = zeros(1,length(xyzdir));
    for i=1:length(xyzdir)-2
        theta(i) = angle3(xyzdir(:,i), xyzdir(:,i+1));
    end
    thetanorm = theta-mean(theta);
    branchdata.Tort_new(b) = std(thetanorm);
       
    J = linspace(0,1,rad_precision+2);
    J = J(2:end-1);
    r = NaN*ones(rad_precision,1);
    a = r;
    e = r;
    o = r;
    k = 1;
    
    for j=J
        % RADIUS
        intp = branchdata.Interp{b};
        tmax = intp.breaks(end);
        basep = ppval(intp,tmax*j);
        nextp = ppval(intp,tmax*j+tmax/pts_per_branch);
        normal = basep-nextp;
        normal = normal/norm(normal);
        
        % Slicing along the normal direction of the vessel
        try
            [sl,xs,ys,zs] = obliqueslice(double(bw),round(basep)',normal',...
                'Method','nearest');
            
            % Looking for the point of the skeleton on the slice
            skp = abs(xs-basep(1))<1 & abs(ys-basep(2))<1 & abs(zs-basep(3))<1;
            
%             NEW APPROACH WITH RECURSION
            [x_skp, y_skp] = ind2sub(size(skp),find(skp));
            truesec = floodimg(sl,[round(mean(x_skp)),round(mean(y_skp))])==2;
            
%             OLD APPROACH WITH GEODESICS
%             geo = bwdistgeodesic(logical(sl),skp);
%             truesec = (~isinf(geo) & ~isnan(geo));

            area = nnz(truesec);
            r(k) = sqrt(area/pi);
            if r(k)==0
                r(k) = 1;
            end
        catch 
            errors = errors+1;
        end
            % LATERAL SURFACE AREA
            a(k)= 2*pi*r(k)*branchdata.Len(b)/numel(J);
            
            % ECCENTRICITY and ORIENTATION
            currecc = regionprops(truesec,'MajorAxisLength','MinorAxisLength','Orientation');
            if ~isempty(currecc)
                e(k) = sqrt(1-currecc.MinorAxisLength^2/currecc.MajorAxisLength^2);
                o(k) = currecc.Orientation;
            end
        k = k+1;
    end
    branchdata.Rad(b) = mean(r,'omitnan');
    branchdata.Alat(b) = mean(a,'omitnan');
    branchdata.Eccent(b) = mean(e,'omitnan');
    branchdata.Orientat(b) = mean(o,'omitnan');
    
    % Approximating the vessel hydraulic resistance
    mu = 1e-3; % Viscosity
    branchdata.Res(b) = 8*mu*branchdata.Len(b)/(pi*branchdata.Rad(b)^4);
    
    % Goodness condition of the vessel: Length > Radius*3 (arbitrary)
    branchdata.isGood(b) = branchdata.Len(b) > branchdata.Rad(b)*3;
    
    fprintf(repmat('\b',1,displine))
    displine = fprintf(strcat(string(b),"/",string(tot_branches)," branches analyzed"));
    
end

fprintf(repmat('\b',1,displine))
if errors>0
   disp(strcat(string(errors),...
       " vessels had to be sliced with less precision than specified in 'radius_precision'")); 
end

% Scaling the metrics based on the downsampling and conversion to
% micrometers
branchdata.Rad = branchdata.Rad*pxdens(1)*downfactor;
branchdata.Len = branchdata.Len*pxdens(1)*downfactor;
branchdata.Alat = branchdata.Alat*(pxdens(1)^2)*(downfactor^2);

% CALCULATING THE REAL LATERAL AREA
realAlat = nnz(bwmorph3(bw,'remove'));

try
    mvn.info.chrono = [mvn.info.chrono; {"Morphological parameters calculation",toc}];
catch
end

clear a area b basep c currecc d e geo intp J k len mu nextp normal r skp...
    sl theta thetanorm truesec v1 xyzpath xyzdir xs ys zs x_skp y_skp i j n...
    tmax displine

%% PCR ANALYSIS
tic
disp("> PCR analysis ");
[x_bw, y_bw, z_bw] = ind2sub([h,w,l],find(bw));
xyzbw = [x_bw y_bw z_bw];
xyzbw_norm = xyzbw - repmat(mean(xyzbw),numel(x_bw),1);
covmat = cov(xyzbw_norm);
[~, eigval] = eig(covmat);
eigval = flipud(diag(eigval));

mvn.pcr = cumsum(eigval)/sum(eigval);

try
    mvn.info.chrono = [mvn.info.chrono; {"PCR analysis",toc}];
catch
end

clear x_bw y_bw z_bw xyzbw xyzbw_nowm covmat eigval

%% .pts CONVERSION FOR CFD INPUT
%-------------------------------------------------------------------------%
% Luca Possenti, Giustina Casagrande, Simone Di Gregorio, Paolo Zunino, 
% Maria Laura Costantino,
% "Numerical simulations of the microvascular fluid balance with a 
% non-linear model of the lymphatic system"
% Microvascular Research, Volume 122, 2019, ISSN 0026-2862,
% https://doi.org/10.1016/j.mvr.2018.11.003.
%-------------------------------------------------------------------------%
if print_pts
    disp("> Creating the input files for the CFD");
    tic
    conversionpath = strcat(string(pathtoimg),"_bifurcation.pts");
    FID = fopen(conversionpath, 'w');
    fprintf(FID,'BEGIN_LIST\n');
    
    s = branchdata.From;
    e = branchdata.To;

    for i = 1 : tot_branches

        fprintf (FID, 'BEGIN_ARC\n');
        fprintf (FID, char(strcat("BC ",string(branchdata.CatFr(i)),"\n")));
        fprintf (FID, char(strcat("BC ",string(branchdata.CatTo(i)),"\n")));
        % Coordinates must be normalized
        fprintf(FID,'%d\t%f\t%f\t%f\tstart\n', i, s(i,1)/max([w,h,l]), ...
            s(i,2)/max([w,h,l]), s(i,3)/max([w,h,l]));
        fprintf(FID,'%d\t%f\t%f\t%f\tend\n', i, e(i,1)/max([w,h,l]), ...
            e(i,2)/max([w,h,l]), e(i,3)/max([w,h,l]));

        for j = 1:pts_per_branch-2
            fprintf(FID,'%d\t%f\t%f\t%f\tpoint\n', i, branchdata.xInt(i,j)/max([w,h,l]),...
                branchdata.yInt(i,j)/max([w,h,l]), branchdata.zInt(i,j)/max([w,h,l]));
        end
        fprintf (FID, 'END_ARC\n');
    end

    fprintf(FID,'END_LIST\n');
    fclose(FID);
    
    conversionpath = strcat(string(pathtoimg),"_radius.pts");
    FID = fopen(conversionpath, 'w');
    r = branchdata.Rad;
    fprintf(FID,'BEGIN_LIST\n');
    for i= 1:tot_branches
        fprintf(FID,'%f\n', r(i)/max([w,h,l])/mean(pxdens(1)/pxdens(2))*downfactor);
    end
    fprintf(FID,'END_LIST\n');
    fclose(FID);

    [ptscheck.iscorrect, ptscheck.where_error] = ...
        checkPts(strcat(string(pathtoimg),"_bifurcation.pts"));
    mvn.info.ptsOK = ptscheck.iscorrect;
    disp(strcat("    PTS check: ",string(ptscheck.iscorrect)));
    if ~ptscheck.iscorrect
        disp(['Error in .pts: ' ptscheck.where_error]);
    end
    try
        mvn.info.chrono = [mvn.inf\o.chrono; {"Generation of .pts file",toc}];
    catch
    end
    
    clear e s end_n start_n FID flag_bpEN flg bpST int p j k xbp ybp zbp xep ...
        yep zep x_tr y_tr z_tr r i appPoint tmax interv ans
end

%% SAVING RESULTS
% Creating the struct 'mvn'
mvn.branchdata = branchdata;
mvn.mRad = mean(branchdata.Rad(:),'omitnan');
bwd = bwdist(~bw);
bwdsk = bwd; bwdsk(sk==0) = 0;
mvn.mRad_REAVER = mean(bwdsk(bwdsk>0),'all')*pxdens(1);
mvn.mLen = mean(branchdata.Len(:),'omitnan');
mvn.mTort = mean(branchdata.Tort(:),'omitnan');
mvn.mEcc = mean(branchdata.Eccent(:),'omitnan');
mvn.volFrac = nnz(bw)/(numel(bw));
mvn.approxAlat = sum(branchdata.Alat(:),'omitnan');
mvn.realAlat = realAlat*pxdens(1)^2*downfactor^2;
mvn.S_over_V = mvn.realAlat/(numel(bw)*pxdens(1)^3*downfactor^3);
mvn.numSubN = numsn;
for i = 1:numsn
    if nnz(mvn.skel.graph.Nodes.subN == i) == 1
        mvn.numSubN = mvn.numSubN - 1;
    end
end

% Completing the struct 'info'
mvn.info.name = pathtoimg;
mvn.info.pxdensity = pxdens;
mvn.info.downfactor = downfactor;
mvn.info.pts_per_branch = pts_per_branch;
mvn.info.ragPrecision = rad_precision;
mvn.info.lungThr = lung_thr;
mvn.info.chrono.Parziale = cumsum(mvn.info.chrono.Duration);
resultsave_path = strcat(string(pathtoimg),append_to_path,".mat");
save(resultsave_path,'mvn');

disp(strcat("> Results saved at: ",resultsave_path));

clear resultsave_path i bwd bwdsk
end
