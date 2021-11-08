function savefigs(varargin)
% SAVEFIGS saves all the figures currently open on the desktop. Figures
% will be saved as JPG image files by default, with a filename
% corresponding to the figure property "Name" settable in MATLAB. The
% extension can be specified as a parameter to the function
%
% Example:
% savefigs;  --> Saves all open figures as JPG
%
% savefigs .png;    --> Saves all open figures as PNG

FolderName = "C:\Users\alber\Desktop\"; % Change this to your desktop path
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
if nargin == 0
    extension = ".jpg";
else
    extension = varargin{1};
end
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = string(get(FigHandle, 'Number'))+"_"+string(get(FigHandle, 'Name'));
    saveas(FigHandle, strcat(FolderName, string(FigName), extension));
end
end