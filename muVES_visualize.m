% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota
%
% For informations please contact:
%   alberto2.rota@mail.polimi.it
%   alberto_rota@outlook.com
% ========================================================================%
%% VISUALIZATION TOOL FOR MICROVASCULAR NETWORKS ANALYZED WITH "NOME"
% Press "Run" to choose which network to visualize. 
[name,path] = uigetfile("*.mat");
load(strcat(path,name));
%==========================================================================%
% Set to 1 the visualizations that you want to be shown, 0 otherwise
flat_img              = 0;  
slices                = 0;  % --> Only for 3D data
segmentation          = 1;  % --> Segmentation and Skeletonization are on 
skeleton              = 1;  %     the same image for 2D data, set either to
interpolation         = 0;  %     1 and both will bw shown
graph                 = 1;
histograms            = 0;
radius                = 1;
lengthh               = 1;
tortuosity            = 1;
pts_classification    = 0;
%==========================================================================%
if numel(size(mvn.bw))==3
    mVN_3D_visualize()
else
    mVN_2D_visualize()
end