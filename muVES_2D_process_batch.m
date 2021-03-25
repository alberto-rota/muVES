% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto2.rota@mail.polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@mail.polimi.it
% ========================================================================%
%% ANALYZE IMAGE BATCH - 2D
if ispc
    folderpath = uigetdir;
    files = dir(folderpath);
    for k=3:length(files)
        pathtoimg=strcat(files(k).folder,"\",files(k).name);
        try
            muVES_2D(pathtoimg);
            splitpath = split(pathtoimg,".");
            load(strcat(splitpath(1),".mat"));
            writetable(mvn.branchdata(:,[1 4 6 13 14 15 16 17]),...
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
            muVES_2D(pathtoimg);
            splitpath = split(pathtoimg,".");
            load(strcat(splitpath(1),".mat"));
            writetable(mvn.branchdata(:,[1 4 6 13 14 15 16 17]),...
                strcat(folderpath,"/batch.xlsx"),'Sheet',string(k-2));
        catch ex
        end
    end
else
    error("Error determining your OS");
end