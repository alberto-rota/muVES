% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota
%
% For informations please contact:
%   alberto2.rota@mail.polimi.it
%   alberto_rota@outlook.com
% ========================================================================%
%% ANALYZE IMAGE BATCH - 3D
if ispc
    folderpath = uigetdir;
    files = dir(folderpath);
    for k=3:length(files)
        pathtoimg=strcat(files(k).folder,"\",files(k).name);
        try
            muVES_3D(pathtoimg);
            splitpath = split(pathtoimg,".");
            load(strcat(splitpath(1),".mat"));
            writetable(mvn.branchdata(:,[1 5 7 14 16 17 18 19]),...
                strcat(folderpath,"\batch.xlsx"),'Sheet',string(k-2));
        catch ex
        end
    end
    
elseif ismac
    folderpath = uigetdir;
    files = dir(folderpath);
    for k=3:length(files)
        pathtoimg=strcat(files(k).folder,"/",files(k).name);
        try
            muVES_3D(pathtoimg);
            splitpath = split(pathtoimg,".");
            load(strcat(splitpath(1),".mat"));
            writetable(mvn.branchdata(:,[1 5 7 14 16 17 18 19]),...
                strcat(folderpath,"/batch.xlsx"),'Sheet',string(k-2));
        catch ex
        end
    end
else
    error("Error determining your OS");
end