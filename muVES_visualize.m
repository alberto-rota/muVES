% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto2.rota@mail.polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@mail.polimi.it
% ========================================================================%
%% VISUALIZATION TOOL FOR MICROVASCULAR NETWORKS ANALYZED WITH muVES
% Press "Run" to choose which network to visualize. 
[name,path] = uigetfile("*.mat");
load(strcat(path,name));
%==========================================================================%
% Set to 1 the visualizations that you want to be shown, 0 otherwise
flat_img              = 1;  
slices                = 1;  % --> Only for 3D data
segmentation          = 1;  % --> Segmentation and Skeletonization are on 
skeleton              = 1;  %     the same image for 2D data, set either to
interpolation         = 1;  %     1 and both will bw shown
graph                 = 1;
histograms            = 1;
radius                = 1;
lengthh               = 1;
tortuosity            = 1;
pts_classification    = 1;
%==========================================================================%
if numel(size(mvn.bw))==3
    muVES_3D_visualize()
else
    muVES_2D_visualize()
end