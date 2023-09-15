% =============================== µVES ================================== %
% Copyrights © 2021     Luca Possenti
%
% For informations please contact:
%   luca.possenti@polimi.it
% 
% ========================================================================%

%% Computing velocity in the networks
% This script compute the velocity of the fluid within the network. The
% network graph is required (see muVES*D.m). 
% BCs are defined by pressure inlet and p = 0 as outlet. Flow rate to be
% implemented. Define BC directly in the code below.
% [pressure] = Pa, [flowrate] = um^3/s, [velocity] = um/s

%% Reading the network
% Press "Run" to choose the network. 
[name,path] = uigetfile("*.mat");
load(strcat(path,name));

%% Viscosity
mu = 1e-3; %[Pa * s]

%% Retrieving BC (p = pressure, q = flow rate (totQ), s = symmetry, o = outlet)
% zone1.type = 'p';
% zone1.value = 50;
% zone2.type = 's';
% zone2.value = 0;
% zone3.type = 'o';
% zone3.value = 50;
% zone4.type = 's';
% zone4.value = 0;

zone1.type = 'p';
zone1.value = 50;
zone2.type = 's';
zone2.value = 5;
zone3.type = 'o';
zone3.value = 50;
zone4.type = 's';
zone4.value = 0;

BC = [zone1 zone2 zone3 zone4];


%% Retrieving branch and nodes info
% Branches are identified by the uniqueID mvn.branchdata.Num
% Nodes should be identified to build a connectivity matrix
uniqueNodes = mvn.branchdata.From(1,:);
for n=2:size(mvn.branchdata.From,1)
    if isempty(find(sum(uniqueNodes==mvn.branchdata.From(n,:),2)==3))
        uniqueNodes = [uniqueNodes; mvn.branchdata.From(n,:)];
    end
end
for n=1:size(mvn.branchdata.To,1)
    if isempty(find(sum(uniqueNodes==mvn.branchdata.To(n,:),2)==3))
        uniqueNodes = [uniqueNodes; mvn.branchdata.To(n,:)];
    end
end
nodeIDs = 1:1:size(uniqueNodes,1);
uniqueNodes = [nodeIDs' uniqueNodes];

%% Building connectivity matrix
connectivity = sparse(size(mvn.branchdata,1),size(uniqueNodes,1)); %branch on the row, nodes on the column
% For each branch find the nodes From (-1) and To (+1)
for b=1:size(mvn.branchdata,1)
    indexFrom = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.From(b,:),'rows'));
    indexTo = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.To(b,:),'rows'));
    connectivity(b,indexFrom) = 1;
    connectivity(b,indexTo) = -1;
end

%% Building resistances matrix [Pa / (um^3/s)]
totalVessels = size(mvn.branchdata,1);
resistance = sparse(totalVessels);
for b=1:totalVessels
    if mvn.branchdata.Rad(b) ==0
        warning('A zero radius has been detected in the vascular network. Radius has been substitued with 1 um.')
        resistance(b,b) = 8 * mu * mvn.branchdata.Len(b) / (pi * (1)^4);
    else
        resistance(b,b) = 8 * mu * mvn.branchdata.Len(b) / (pi * mvn.branchdata.Rad(b)^4);
    end
end

%% Plotting
figure 
hold on
for b=1:totalVessels
    %indexFrom = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.From(b,:),'rows'));
    %indexTo = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.To(b,:),'rows'));
    coord = [mvn.branchdata.From(b,:); mvn.branchdata.To(b,:)];
    plot3(coord(:,1),coord(:,2),coord(:,3));
    if contains(char(mvn.branchdata.CatFr(b)),'DIR')
        plot3(coord(1,1),coord(1,2),coord(1,3),'ro');
    elseif contains(char(mvn.branchdata.CatFr(b)),'INT')
        plot3(coord(1,1),coord(1,2),coord(1,3),'bo');
    else
        plot3(coord(1,1),coord(1,2),coord(1,3),'ko');
    end
    if contains(char(mvn.branchdata.CatTo(b)),'DIR')
        plot3(coord(2,1),coord(2,2),coord(2,3),'ro');
    elseif contains(char(mvn.branchdata.CatTo(b)),'INT')
        plot3(coord(2,1),coord(2,2),coord(2,3),'bo');
    else
        plot3(coord(2,1),coord(2,2),coord(2,3),'ko');
    end
end
title('Network geometry (nodes: red-pressureBC, blue-junctions, black-noFlowBC)')
hold off

%% Building the system
A = [connectivity resistance; sparse(size(uniqueNodes,1),size(uniqueNodes,1)) connectivity'];
B = zeros(size(mvn.branchdata,1)+size(uniqueNodes,1),1);

%% Implementing BC
for b=1:totalVessels
    % nodes from
    indexFrom = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.From(b,:),'rows'));
    if contains(char(mvn.branchdata.CatFr(b)),'DIR')
        zone = strsplit(char(mvn.branchdata.CatFr(b)));
        zone = str2double(zone{2});  
        if BC(zone).type == 'p'
            A(totalVessels+indexFrom,:)= 0;
            A(totalVessels+indexFrom,indexFrom) = 1;
            B(totalVessels+indexFrom) = BC(zone).value;
        elseif BC(zone).type == 'o'
            A(totalVessels+indexFrom,:)= 0;
            A(totalVessels+indexFrom,indexFrom) = 1;
            B(totalVessels+indexFrom) = 0;
        elseif BC(zone).type == 's'
            % nothing to be done (Q = 0)
        else
            str = strcat('Unknown boundary condition: ',BC(zone).type);
            err(str)
        end
    end
    % nodes to
    indexTo = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.To(b,:),'rows'));
    if contains(char(mvn.branchdata.CatTo(b)),'DIR')
        zone = strsplit(char(mvn.branchdata.CatTo(b)));
        zone = str2double(zone{2});  
        if BC(zone).type == 'p'
            A(totalVessels+indexTo,:)= 0;
            A(totalVessels+indexTo,indexTo) = 1;
            B(totalVessels+indexTo) = BC(zone).value;
        elseif BC(zone).type == 'o'
            A(totalVessels+indexTo,:)= 0;
            A(totalVessels+indexTo,indexTo) = 1;
            B(totalVessels+indexTo) = 0;
        elseif BC(zone).type == 's'
            % nothing to be done (Q = 0)
        else
            str = strcat('Unknown boundary condition: ',BC(zone).type);
            err(str)
        end
    end
end 

%% Solving
X = A\B;
flowrate=X(end-totalVessels+1:end); %um^3/s
flowrate_ul = flowrate*1e-9; %ul/s
pressure=X(1:end-totalVessels); %Pa

%% Computing velocity
radius = mvn.branchdata.Rad;
if nnz(radius)~=length(radius)
    %warning already done before
    mask = radius==0;
    radius(mask) = 1;
end
velocities = flowrate./(pi.*radius.^2); %um/s

%% Plotting
% Flowrate
minQ = min(abs(flowrate_ul));
maxQ = max(abs(flowrate_ul));
if minQ == maxQ
    minQ = minQ*0.9;
    maxQ = maxQ*1.1;
end
figure
title('Flow rate (\mul/s)')
c=colormap(jet);
hold on
for b=1:totalVessels
    coord = [mvn.branchdata.From(b,:); mvn.branchdata.To(b,:)];
    col = ceil((abs(flowrate_ul(b))-minQ)/(maxQ-minQ)*255)+1;
    if isnan(col)
        warning('Color not defined. Check data carefully!')
        continue;
    end
    plot3(coord(:,1),coord(:,2),coord(:,3),'Color',c(col,:),'LineWidth',7);
end
hold off
colorbar
caxis([minQ maxQ])

% velocity
minV = min(abs(velocities));
maxV = max(abs(velocities));
if minV == maxV
    minV = minV*0.9;
    maxV = maxV*1.1;
end
figure
title('Velocity (\mum/s)')
c=colormap(jet);
hold on
for b=1:totalVessels
    coord = [mvn.branchdata.From(b,:); mvn.branchdata.To(b,:)];
    col = ceil((abs(velocities(b))-minV)/(maxV-minV)*255)+1;
    if isnan(col)
        warning('Color not defined. Check data carefully!')
        continue;
    end
    plot3(coord(:,1),coord(:,2),coord(:,3),'Color',c(col,:),'LineWidth',7);
end
hold off
colorbar
caxis([minV maxV])

% Pressure
minP = min(abs(X(1:end-totalVessels)));
maxP = max(abs(X(1:end-totalVessels)));
if minP == maxP
    minP = minP*0.9;
    maxP = maxP*1.1;
end
figure
title('Pressure (Pa)')
c=colormap(jet);
hold on
for b=1:totalVessels
    coord = [mvn.branchdata.From(b,:); mvn.branchdata.To(b,:)];
    indexFrom = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.From(b,:),'rows'));
    indexTo = find(ismember(uniqueNodes(:,2:4),mvn.branchdata.To(b,:),'rows'));
    if isnan(indexFrom) || isnan(indexTo)
        continue;
    end
    colFrom = ceil((abs(X(indexFrom))-minP)/(maxP-minP)*255)+1;
    colTo = ceil((abs(X(indexTo))-minP)/(maxP-minP)*255)+1;
    if isnan(colFrom) || isnan(colTo)
        warning('Color not defined. Check data carefully!')
        continue;
    end
    plot3(coord(:,1),coord(:,2),coord(:,3),'k');
    plot3(coord(1,1),coord(1,2),coord(1,3),'-o','Color','k','MarkerSize',10,'MarkerFaceColor',c(colFrom,:));
    plot3(coord(2,1),coord(2,2),coord(2,3),'-o','Color','k','MarkerSize',10,'MarkerFaceColor',c(colTo,:));
end
hold off
colorbar
caxis([minP maxP])

%% Saving
str = strcat(path,'pressureAndVelocity_',name);
save(str,'flowrate','flowrate_ul','pressure','velocities')