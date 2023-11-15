% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto1.rota@polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@polimi.it
% ========================================================================%
% SPECIFY THE SETTINGS IN THE "muVES settings.txt" FILE 
% ======================================================================= %
%% SETUP AND IMAGE LOADING
function mvn = muVES_2D(varargin)
% ======================================================================= %
% REQUIREMENTS: 
req = {'Image Processing Toolbox', 'Curve Fitting Toolbox'};
% ======================================================================= %
if nargin == 0
% Specify the path to the file. The extension must be separately specified
% in 'extension'. If those fields are left empty ('') or if the path is not
% correct, a window for file selection will be opened.
% Example:      pathtoimg = 'C:\...\...\myfolder\myfile';
%               extension = '.oib'; 
pathtoimg = "";
extension = ".tif";

    try
        img = imread(strcat(pathtoimg,extension));
    catch
        [name, folder] = uigetfile("*.*");
        path = strcat(folder,name);
        splpath = strsplit(path,".");
        pathtoimg = splpath(1);
        extension = strcat(".",splpath(2));
    end
else
    splpath = strsplit(varargin{1},".");
    pathtoimg = splpath(1);
    extension = strcat(".",splpath(2));
end

fid = fopen("muVES settings.txt");
% =========================================================================%
% All the results are stored in the struct 'mvn', accessible in 
% the 'filename.mat' created in the same folder as the input.
% --> See 'MVN legend.txt' for all the details of the values in 'mvn'. 
% ========================================================================%
% 'pxdens' contains the space resolution of the microscope in the 3
% dimensions, specified in micrometers/pixel. The algorithm interpolates
% the volume so that the resolution in the Z direction is the same to the
% resoution in the X and Y direction
xds = strsplit(fgetl(fid),":");
yds = strsplit(fgetl(fid),":");
zds = strsplit(fgetl(fid),":");
pxdens = [str2double(xds{2}) str2double(yds{2})];
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
rag_precision = str2double(pds{2});   % [POSITIVE INTEGER SCALAR]
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
tic
switch extension
    case {'.oib', '.n2d', '.tif'}
        pathtoimg = char(pathtoimg);
        disp(strcat("Processing: ",string(pathtoimg)));
%         try
%-------------------------------------------------------------------------%
% Melissa Linkert, Curtis T. Rueden, Chris Allan, Jean-Marie Burel, 
% Will Moore, Andrew Patterson, Brian Loranger, Josh Moore, Carlos Neves,
% Donald MacDonald, Aleksandra Tarkowska, Caitlin Sticco, Emma Hill, Mike
% Rossner, Kevin W. Eliceiri, and Jason R. Swedlow (2010) 
% "Metadata matters: access to image data in the real world." 
% The Journal of Cell Biology 189(5), 777-782. doi: 10.1083/jcb.201004104
%-------------------------------------------------------------------------%
            img_3d = bfOpen3DVolume(char(strcat(pathtoimg,extension)));
            img = img_3d{1,1}{1,1};
            splpath = strsplit(img_3d{1}{2},".");
            pathtoimg = splpath(1);
            extension = strcat(".",splpath(2));
            img = mat2gray(max(img, [], 3));
            mvn.info.hashtable = string(img_3d{1,2});
%         catch
%             error('Path Invalido o il file non è nel formato corretto');
%         end
    case {'.bmp'}
        img = imread(strcat(pathtoimg,extension));
        img = img(:,:,1);
    otherwise
        warning(strcat("Unsupported file type: ",pathtoimg,extension));
        return;
end
tic
disp(strcat("> Processing: ",string(pathtoimg),extension));

% Inizializzazione della table con i tempi
mvn.info.chrono = table("File reading and initialization",toc,'VariableNames',...
    {'Phase','Duration'});
mvn.raw = img;

clear path splpath xds yds zds pds fid

% DOWNAMPLING del fattore 'downfactor' specificato (per comodità dichiarato
% nella sezione 1). Il downsampling viene effettuato perchè MATLAB
% richiede tempi -relativamente- lunghi per ricavare
% isosurfaces da una grande quantità di dati. Con downfactor = 2-3 il
% risultato è comunque molto nitido. Al posto di considerare tutti i voxel
% nelle 3 dimensioni, vengono "campionati" con una distanza specificatela dal
% valore di 'downfactor'
img = img(1:downfactor:end, 1:downfactor:end, 1:downfactor:end);

%% IMAGE ENHANCEMENT AND BINARIZATION
tic
disp("> Enhancing the image and performing binarization ");
img_mf = medfilt2(img, [10 10]);  % Median filtering
gamma = 0.5;
img_en = imadjust(img_mf, [0;1], [0;1], gamma); % Dark pixel contrast enhancement
bw = imbinarize(img_en); %bw = medfilt2(bw, [5 5]);% Removes salt&pepper 
bw = bwmorph(bw,'spur');
bw = bwmorph(bw,'clean')
bw = imgaussfilt(double(bw),5);
bw = imbinarize(bw,0.5);
mvn.bw = bw;

try
    mvn.info.chrono = [mvn.info.chrono; {"Image Enhancement e Segmentazione",toc}];
catch
end

clear ca

%% SKELETONIZATION 
% La funzione Skeleton3D richiede un'immagine già segmentata, quindi una
% matrice 3D di valori logici (bw)
tic
disp("> Computing the skeleton ");
% 'sk' ha le stesse dimensioni della matrice di partenza
sk = bwskel(bw)';
% Le nuove dimensioni calcolate sono dipendenti dal valore di 'downfactor'
h=size(sk,1);
w=size(sk,2);

mvn.skel.sk = sk;

try
    mvn.info.chrono = [mvn.info.chrono; {"Scheletrizzazione con bwmorph",toc}];
catch
end

clear ca

%% COMPUTING BRANCHPOINTS AND ENDPOINTS
% Le funzioni bwmorph3 per endpoints e branchpoints sono funzioni
% precodificate in MATLAB.
tic
disp("> Finding branchpoints ");
% Come prima, 'bpoints' è una matrice di valori logici della stessa
% dimensione dell'imagine principale (1 dove c'è il branchpoint, zero dove
% non c'è il branchpoint. Le coordinate di questi valori logici
% devono essere riportate nei vettori 'x_bp', 'y_bp' e 'z_bp'. Discorso
% analogo vale per 'epoints'
bpoints = bwmorph(sk,'branchpoints');
epoints = bwmorph(sk,'endpoints');
mvn.skel.bp = bpoints;
mvn.skel.ep = epoints;

try
    mvn.info.chrono = [mvn.info.chrono; {"Calcolo di Endpoints e Branchpoints",toc}];
catch
end

%% INTERPOLATION AND NODE CLASSIFICATION
tic
disp("> Interpolating the branches ");
% Suddivisione in branches con la funzione Skel2Graph3D. Le informazioni
% sui vari branches sono contenuti nel vettore di struct 'link'. Le
% informazioni sui vari nodi sono contenuti nel vettore di struct 'node'.
% L'interpolazione ci permette di corregge il problema dei branchpoint
% tripli
%-------------------------------------------------------------------------%
% Copyright (c) 2016, Philip Kollmannsberger
% All rights reserved.

% Kerschnitzki, Kollmannsberger et al.,
% "Architecture of the osteocyte network correlates with bone material 
% quality."
% Journal of Bone and Mineral Research, 28(8):1837-1845, 2013.
%-------------------------------------------------------------------------%
[adj, node, link] = Skel2Graph3D(sk,lung_thr);
    
    % Il numero di branches è dato dalla quantià di elementi di 'link'
    tot_branches = numel(link);
    
    % Viene creata una table 'vessel data' che contiene le seguenti
    % informazioni:
    %   -'Number' è il numero identificativo del vaso in questione
    %   -'xPath' contiene le coordinate X dei punti che appartengono al branch
    %   -'yPath' contiene le coordinate Y dei punti che appartengono al branch
    %   -'zPath' contiene le coordinate Z dei punti che appartengono al branch
    %   -'Interp' contiene i parametri (plottabili con 'fnplt') dell'interpolazione
    branchdata = table(0,{0},{0},[0,0],{0},[0,0],{0},{0},...
        zeros(1,pts_per_branch-2),zeros(1,pts_per_branch-2),...
        'VariableNames', {'Num','xPath','yPath','From','CatFr','To','CatTo',...
        'Interp','xInt','yInt'});
    
    % Di ogni branch vengono estratte le coordinate in 'path_x','path_y' e
    % 'path_z', aggiunte alla table e poi plottate. Al grafico vengono poi
    % aggiunti i nodi dei vari branches, anch'essi savati nella table.
    for b=1:tot_branches
        path_idx = link(b).point;
        [path_x, path_y] = ind2sub([h,w],path_idx);
        
        interp = cscvn([path_x; path_y]);
        x_from = node(link(b).n1).comx;
        y_from = node(link(b).n1).comy;
        x_to = node(link(b).n2).comx;
        y_to = node(link(b).n2).comy;
        from = [x_from,y_from];
        to = [x_to,y_to];
        toAdd = {b, path_x, path_y, from,"", to, "",interp, ...
            zeros(1,pts_per_branch-2),zeros(1,pts_per_branch-2)};
        branchdata = cat(1,branchdata,toAdd);
    end
    branchdata(1,:) = [];
    
    clear from to toAdd curr_pp  path_x path_y path_z a c
    
    % CONVERSIONE IN GRAFO E CLASSIFICAZIONE DEI NODI
    % Grazie alla matrice di adiacenza trovata con Skel2Graph3D si ricava la
    % conversione dello scheletro in un grafo. Ogni nodo è identificato dal
    % proprio indice e dalle tre coordinate spaziali
    tic
    G = graph(adj);
    for n = 1:numel(node)
        G.Nodes.x(n) = node(n).comx;
        G.Nodes.y(n) = node(n).comy;
        G.Nodes.subN(n) = 0;
    end
    

    % Con la funzione flood si distinguono diverse sottoreti disconnesse nel
    % grafo. Ad ogni nodo della rete viene assegnato il valore 'subN',
    % identificativo della sottorete di appartenenza. Si veda la descrizione
    % della funzione flood nell'ultima sezione.
    numsn = 1;
    while any(G.Nodes.subN == 0)
        non_lab = find(G.Nodes.subN==0);
        G = floodgraph(G,non_lab(1),numsn);
        numsn=numsn+1;
    end
    mvn.skel.graph = G;
    
 if print_pts
    % Si crea, per ogni branch, un vettore di punti in quantità definita da
    % 'pts_per_branch', uguale per ogni ramo indipendentemente dalla lunghezza
    for i=1:tot_branches
        tmax = max(branchdata.Interp{i}.breaks);
        interv = linspace(tmax/(pts_per_branch-2),tmax-tmax/(pts_per_branch-2),...
            (pts_per_branch-2));
        for j = 1:pts_per_branch-2
            p = ppval(branchdata.Interp{i}, interv);
            branchdata.xInt(i,:) = p(1,:);
            branchdata.yInt(i,:) = p(2,:);
        end
    end
    
    % Estrazione delle coordinate dei punti di start e end di ogni branches
    s = branchdata.From;
    e = branchdata.To;
    
    % Ogni nodo dello scheletro viene categorizzato con la condizione opportuna
    % per la simulazione fluidodinamica:
    % -"INT" per i branchpoints
    % -"DIR" per gli endpoints vicino ai bordi, con un valore di pressione
    %  dipendente dal bordo più vicino (da 1 a 4)
    % -"MIX" per gli endpoints non vicino ai bordi
    G.Nodes.deg = degree(G);
    
    for i = 1: tot_branches
        % Individuazione dei branchpoints della rete:
        % Un branchpoint è un nodo di grado 1 nel grafo
        if G.Nodes.deg(all(branchdata.From(i,:) == table2array(G.Nodes(:,1:2)),2)) == 1
            flag_bpST = 0;
        elseif G.Nodes.deg(all(branchdata.From(i,:) == table2array(G.Nodes(:,1:2)),2)) > 1
            flag_bpST = 1;
        end
        
        % Se lo start è un punto di branch ha bisogno di una condizione "INT",
        % de è un punto di end ha bisogno di una condizione "MIX" o "DIR"
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
        
        if G.Nodes.deg(all(branchdata.To(i,:) == table2array(G.Nodes(:,1:2)),2)) == 1
            flag_bpEN = 0;
        elseif G.Nodes.deg(all(branchdata.To(i,:) == table2array(G.Nodes(:,1:2)),2)) > 1
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
        idx = find(all(branchdata.From(i,:) == [G.Nodes.x G.Nodes.y],2));
        branchdata.subN(i) = G.Nodes.subN(idx);
    end
    
    % Il primo vaso nel file .pts deve contenere una condizione "DIR".
    % Viene quindi effettuato uno scambio nell'ordine dei vasi
    for i=1:tot_branches
        if branchdata.CatFr(i) == "DIR 1" || branchdata.CatTo(i) == "DIR 1" || ...
                branchdata.CatFr(i) == "DIR 2" || branchdata.CatTo(i) == "DIR 2" ||...
                branchdata.CatFr(i) == "DIR 3" || branchdata.CatTo(i) == "DIR 3" ||...
                branchdata.CatFr(i) == "DIR 4" || branchdata.CatTo(i) == "DIR 4"
            branchdata([1,i],:) = branchdata([i,1],:);
            break;
        end
    end
    
    % Per il ptimo vaso si mette una condizione "DIR" sul punto di start. Per
    % tutti gli altri vasi, la condizione "DIR" deve essere sull'end. Viene
    % quindi invertito l'ordine dei punti dello scheletro (che devono
    % necessariamente andare da start a end e non viceversa)
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
    end
    
    % Per il resto dei vasi nella rete (dal secondo all'ultimo), la condizione
    % "DIR" (se presente) deve essere sull'end. Viene effettuata la stessa
    % operazione di scambio
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
            branchdata(branchdata.CatFr == "MIX" & cond,:) = [];
        end
    end
% %     tot_branches = size(branchdata,1);
 end
 sk = zeros(h,w);
 for i=1:size(branchdata,1)
     for j=1:numel(branchdata.xPath{i})
         sk(branchdata.xPath{i}(j),branchdata.yPath{i}(j)) = 1;
     end
 end
 bpoints = bwmorph(sk,'branchpoints');
 epoints = bwmorph(sk,'endpoints');
 [x_sk, y_sk] = ind2sub([h,w],find(sk));
 [x_bp, y_bp] = ind2sub([h,w],find(bpoints));
 [x_ep, y_ep] = ind2sub([h,w],find(epoints));
 
 branchdata.zInt = zeros(size(branchdata,1));
 branchdata.From = [branchdata.From zeros(size(branchdata,1),1)];
 branchdata.To = [branchdata.To zeros(size(branchdata,1),1)];

try
    mvn.info.chrono = [mvn.info.chrono; {"Interpolazione e refining",toc}];
catch
end

clear non_lab ph x_from y_from z_from x_to y_to z_to

%% MORPHOLOGICAL MEASUREMENTS
tic
disp("> Morphological Measurements ");
displine = 0;
errors = 0;
% Per il calcolo della lunghezza dei branches, si lavora sullo scheletro
% discreto. In particolare, si valuta il numero di voxel adiacenti in
% orizzontale o verticale e il numero di voxel adiacenti in diagonale. La
% lunghezza in verticale/orizzontale è pari a 1 pixel, in diagonale planare
% è pari a 1.41 pixel, in diagonale 3D è pari a 1.73 pixel.
% figure; imagesc(bw); hold on; colormap gray;axis equal;
for b=1:tot_branches
    % Di ogni voxel, viene calcolata la distanza con il successivo. Se la
    % distanza (al quadrato, per evitare di usare l'operazione di radice,
    % richiede molta potenza di calcolo. Non ci cambia niente) è pari a 1
    % l'adiacenza è orizzontale/vericale, altrimenti diagonale
    
    xypath = [branchdata.xPath{b};branchdata.yPath{b}];
    branchdata.Len(b) = sum(sqrt(sum(diff(xypath,1,2).^2)));
    
    % La tortuosità è definita come la lunghezza del vaso (precedentement
    % calcolata) divisa per la distanza tra i punti di inizio e fine.
    % Poichè alcuni la posizione di alcuni nodi viene interpolata (problema
    % dei tripli nodi), può succedere che il valore della lunghezza sia
    % leggermente inferiore alla distanza inizio-fine, risultando in un
    % valore di tortuosità inferiore a 1 (impossibile). Poichè ciò si
    % verifica quasi esclusivamente per vasi corti e diagonali, è lecito
    % costringere il valore della tortuosità sopra 1.

    d = sqrt(sum((branchdata.From(b,:)-branchdata.To(b,:)).^2));
    branchdata.Tort(b) = constrain(branchdata.Len(b)/d,[1 Inf]);
    
    % Modificando il valore di si può valutare il raggio con più precisione
    J = linspace(0,1,rag_precision+2);
    J = J(2:end-1);
    r = NaN*ones(rag_precision,1);
    a = r;
    k = 1;
    
    for j=J
        % Il raggio viene calcolato con una slice secondo la direzione normale
        % alla tangente del branch in un punto scelto come base. La
        % procedura viene effettuata 3 volte, per 3 punti diversi sulla
        % lunghezza e poi viene calcolata la media
        intp = branchdata.Interp{b};
        tmax = intp.breaks(end);
        % La normale di taglio viene scelta tra il punto base ed uno vicino di
        % un fattore 1/pts_per_branch
        bp = ppval(intp,tmax*j);
        np = ppval(intp,tmax*j+tmax/pts_per_branch);
        
        m = -(np(1)-bp(1))/(np(2)-bp(2));
        
        % Slicing nel punto medio sulla direzione normale appena calcolata
        if m<-30 || m>30
            slice = double(bw(:,round(bp(1))))';
            truesec = floodimg(slice',[1 floor(bp(2))]);
            r(k) = nnz(truesec == 2)/2;
        else
            q = bp(2)-m*bp(1);
            slicefun = @(x) m*x+q;
            slicefrom = [1, slicefun(1)];
            sliceto = [h, slicefun(h)];
            [xslice, yslice, slice] = ...
                improfile(bw,[slicefrom(1), sliceto(1)], [slicefrom(2) sliceto(2)]);
            xslice(isnan(slice)) = [];
            yslice(isnan(slice)) = [];
            slice(isnan(slice)) = [];
            if abs(m)<0.3
                [~,skp] = min(abs(xslice-bp(1)));
            else
                [~,skp] = min(abs(yslice-bp(2)));
            end
            truesec = floodimg(slice', [1 skp]);
% UN-COMMENT PER VEDERE IN AZIONE IL CALCOL0 DEL RAGGIO
%             figure;
%             imagesc(bw);hold on;fplot(slicefun);
            truesec = truesec == 2;
            if nnz(truesec)>0
             r(k) = sqrt((max(xslice(truesec))-min(xslice(truesec))).^2+...
                (max(yslice(truesec))-min(yslice(truesec))).^2)/2;
            else
                r(k) = 1;
            end
% UN-COMMENT PER VEDERE IN AZIONE IL CALCOL0 DEL RAGGIO
%             fnplt(branchdata.Interp{b},'b',2);
%             scatter(xslice(truesec),yslice(truesec),'.r');
%             text(mean(xslice(truesec)),mean(yslice(truesec)),string(round(r(k))));
%             close;
            if r(k)==0
                r(k) = 1;
            end
        end
        % Area laterale calcolata cilindica come la somma delle aree
        % laterali su tre segmenti di uguale lunghezza ma di raggio diverso
        a(k)= 2*pi*r(k)*branchdata.Len(b)/numel(J);
        k = k+1;
    end
    mu = 1e-3;
    branchdata.Rad(b) = mean(r,'omitnan');
    branchdata.Alat(b) = mean(a,'omitnan');

    branchdata.Res(b) = 8*mu*branchdata.Len(b)/(pi*branchdata.Rad(b)^4);
    branchdata.check(b) = branchdata.Len(b) > branchdata.Rad(b)*3;
    
    fprintf(repmat('\b',1,displine))
    displine = fprintf(strcat(string(b),"/",string(tot_branches)," branches analyzed"));
end
fprintf(repmat('\b',1,displine))
% Scalamento delle metriche in funzione della densità in pixel e del
% downsampling
branchdata.Rad = branchdata.Rad*pxdens(1)*downfactor;
branchdata.Len = branchdata.Len*pxdens(1)*downfactor;
branchdata.Alat = branchdata.Alat*(pxdens(1)^2)*(downfactor^2);

% Calcolo dell'area laterale complessiva eliminando tutti i voxel "interni"
% (meno di 6-connected) e calcolo di S/V
realAlat = nnz(bwmorph(bw,'remove'));

try
    mvn.info.chrono = [mvn.info.chrono; {"Calcolo delle metriche in 3D",toc}];
catch
end

% clear currcolor_idx dist_to_next xs ys zs sl geo truesec area skp normal...
%     basep nextp intp tmax i j k

%

%% .pts CONVERSION FOR CFD INPUT
if print_pts
tic
%-------------------------------------------------------------------------%
% Luca Possenti, Giustina Casagrande, Simone Di Gregorio, Paolo Zunino, 
% Maria Laura Costantino,
% "Numerical simulations of the microvascular fluid balance with a 
% non-linear model of the lymphatic system"
% Microvascular Research, Volume 122, 2019, ISSN 0026-2862,
% https://doi.org/10.1016/j.mvr.2018.11.003.
%-------------------------------------------------------------------------%
disp("> Creating the .pts file ");
% Path dove verranno salvati i files .pts
conversionpath = strcat(string(pathtoimg),"_bifurcation.pts");
% Apertura file
FID = fopen(conversionpath, 'w');
fprintf(FID,'BEGIN_LIST\n');

% Vengono aggiornati 's' ed 'e' con le coordinate dei punti dopo tutti gli
% scambi effettuati precedentemente
s = branchdata.From;
e = branchdata.To;

for i = 1 : tot_branches
    % Vengono considerati nella simulazione solamente rami con lunghezza
    % superiore ad una soglia impostata dall'utente nella prima sezione
    % dello script
    fprintf (FID, 'BEGIN_ARC\n');
    fprintf (FID, char(strcat("BC ",string(branchdata.CatFr(i)),"\n")));
    fprintf (FID, char(strcat("BC ",string(branchdata.CatTo(i)),"\n")));
    % Le coordinate stampate sul file devono essere normalizzate rispetto
    % alla dimensione maggiore dell'immagine (il massimo tra altezza,
    % larghezza e profondità)
    fprintf(FID,'%d\t%f\t%f\t%f\tstart\n', i, s(i,1)/max([w,h]), ...
        s(i,2)/max([w,h]), s(i,3)/max([w,h]));
    fprintf(FID,'%d\t%f\t%f\t%f\tend\n', i, e(i,1)/max([w,h]), ...
        e(i,2)/max([w,h]), e(i,3)/max([w,h]));
    % 'pts_per_branch' è un valore fisso per ogni branch che identifica
    % quanti nodi considerare in ogni branch. Viene settato all'inizio del
    % codice
    for j = 1:pts_per_branch-2
        fprintf(FID,'%d\t%f\t%f\t%f\tpoint\n', i, branchdata.xInt(i,j)/max([w,h]),...
            branchdata.yInt(i,j)/max([w,h]), branchdata.zInt(i,j)/max([w,h]));
    end
    fprintf (FID, 'END_ARC\n');
end

fprintf(FID,'END_LIST\n');
fclose(FID);
% Su un file differente devono essere stampati, in ordine, i raggi
% corrispondenti ad ogni rampo
conversionpath = strcat(string(pathtoimg),"_radius.pts");
FID = fopen(conversionpath, 'w');
r = branchdata.Rad;
fprintf(FID,'BEGIN_LIST\n');
% Anche il valore del raggio deve essere normalizzato della la stessa
% quantità per cui erano stati normalizzate le coordinate dei punti dello
% scheletro
for i= 1:tot_branches
    fprintf(FID,'%f\n', r(i)/max([w,h])/mean(pxdens(1)/pxdens(2))*downfactor);
end
fprintf(FID,'END_LIST\n');
fclose(FID);

[iscorrect, ~] = checkPts(strcat(string(pathtoimg),"_bifurcation.pts"));
mvn.info.ptsOK = iscorrect;
disp(strcat("    PTS check: ",string(iscorrect)));

clear e s end_n start_n FID flag_bpEN flg bpST int p j k xbp ybp zbp xep ...
    yep zep x_tr y_tr z_tr r i appPoint tmax interv

try
    mvn.info.chrono = [mvn.info.chrono; {"Conversione in .pts",toc}];
catch
end
end
%% SAVING RESULTS
disp("> Saving the results ");
% Viene creata una struct 'mvn' (MicroVascular Network), contenente i dati
%   salienti della rete
mvn.branchdata = branchdata;
bwd = bwdist(~bw);
bwdsk = bwd; bwdsk(imrotate(sk,90)==0) = 0;
mvn.mRad_REAVER = mean(bwdsk(bwdsk>0),'all')*pxdens(1);
mvn.mRad = mean(branchdata.Rad(:),'omitnan');
mvn.mLen = mean(branchdata.Len(:),'omitnan');
mvn.mTort = mean(branchdata.Tort(:),'omitnan');
mvn.areaFrac = nnz(bw)/(numel(bw));
mvn.approxAlat = sum(branchdata.Alat(:),'omitnan');
mvn.realAlat = realAlat*pxdens(1)^2*downfactor^2;
mvn.S_over_V = mvn.realAlat/(numel(bw)*pxdens(1)^3*downfactor^3);
% mvn.numSubN = numsn;
% for i = 1:numsn
%     if nnz(mvn.skel.graph.Nodes.subN == i) == 1
%         mvn.numSubN = mvn.numSubN - 1;
%     end
% end

% Viene creata la struct 'info'
mvn.info.name = pathtoimg;
mvn.info.pxdensity = pxdens;
mvn.info.downfactor = downfactor;
mvn.info.pts_per_branch = pts_per_branch;
mvn.info.ragPrecision = rag_precision;
mvn.info.lungThr = lung_thr;
mvn.info.chrono.Parziale = cumsum(mvn.info.chrono.Duration);
resultsave_path = strcat(string(pathtoimg),append_to_path,".mat");
save(resultsave_path,'mvn');
end
