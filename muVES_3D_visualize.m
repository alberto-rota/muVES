
% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto1.rota@polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@mail.polimi.it
% ========================================================================%
%% VISUALIZATION OF 3D MICROVASCULAR NETWORK
tot_branches = size(mvn.branchdata,1);
vol = mvn.raw;
[h,w,l] = size(mvn.bw);
[x_sk,y_sk,z_sk]=ind2sub([h,w,l],find(mvn.skel.sk));
[x_ep,y_ep,z_ep]=ind2sub([h,w,l],find(mvn.skel.bp));
[x_bp,y_bp,z_bp]=ind2sub([h,w,l],find(mvn.skel.ep));

%% 2D PROJECTION
if flat_img
    figure('Name','2D Projection');
    imagesc(mvn.flat);
    ca = gca;
    ca.Title.String = '2D projection';
end
%% SLICE VIEWER
if slices
    sliceViewer(vol);
end
%% SEGMENTATION
if segmentation
    figure('Name','Segmentation');
    vhandle = volshow(double(mvn.bw));
    vhandle.Renderer = 'Isosurface';
    vhandle.IsosurfaceColor = [1 0 0];
    vhandle.BackgroundColor = [0.4 0.4 0.4];
end
%% SKELETON
if skeleton
    figure('Name','3D Discrete skeleton');
    plot3(x_sk,y_sk,z_sk,'square','Markersize',2,'MarkerFaceColor',[1 1 1],...
        'Color',[1 1 1]);
    hold on
    scatter3(x_ep,y_ep,z_ep,'filled','CData',[87 192 255]./255);
    scatter3(x_bp,y_bp,z_bp,'filled','CData',[255 0 0]./255);
    hold off
    
    set(gca,'Color',[0.2 0.2 0.2]);
    ca = gca;
    ca.Title.String = 'Skeleton';
    daspect([1 1 1]);
    view(3);
end
%% INTERPOLATION
if interpolation
    figure('Name','Interpolated skeleton');
    hold on;
    set(gca,'Color',[0.2 0.2 0.2]);
    title('Interpolated Skeleton');
    view(3);
    daspect([1 1 1]);
for i=1:tot_branches
    fnplt(mvn.branchdata.Interp{i},'y',2);
end
from = mvn.branchdata.From; to = mvn.branchdata.To;
scatter3(from(:,1),from(:,2),from(:,3),'filled','CData',[87 192 255]./255);
scatter3(to(:,1),to(:,2),to(:,3),'filled','CData',[255 0 0]./255);
end
%% GRAPH
if graph 
    figure('Name','Graph conversion');
    ph = plot(mvn.skel.graph,'XData',mvn.skel.graph.Nodes.x,'YData',...
        mvn.skel.graph.Nodes.y,'ZData',mvn.skel.graph.Nodes.z);
    title('Graph');
    try
        for sn=1:max(mvn.skel.graph.Nodes.subN)
            highlight(ph,mvn.skel.graph.Nodes.subN == sn,'NodeColor',abs(rand([1 3])),...
                        'MarkerSize',4);
        end
    catch 
    end
end
%% HISTOGRAMS
if histograms
    % ISTOGRAMMA: Lenhezze
    figure('Name','Histgrams: mesurement distribution for different metrics');
    subplot(2,2,1);
    histogram(mvn.branchdata.Len(:),'FaceColor','r','Normalization','probability');
    xline(mean(mvn.branchdata.Len(:)),'k-.','LineWidth',2);
    title('Length');
    xlabel('Length [\mum]');
    ylabel('Occurrences [%]');
    % ISTOGRAMMA: Tortuosità
    subplot(2,2,2);
    histogram(mvn.branchdata.Tort(:),'FaceColor','b','Normalization','probability');
    xline(mean(mvn.branchdata.Tort(:)),'k-.','LineWidth',2);
    xlabel('Tortuosity [adim.]');
    ylabel('Occurrences');
    title('Tortuosity');
    % ISTOGRAMMA: Raggio
    subplot(2,2,3);
    histogram(mvn.branchdata.Rad(~isnan(mvn.branchdata.Rad)),'FaceColor','g',...
        'Normalization','probability');
    xline(mean(mvn.branchdata.Rad(~isnan(mvn.branchdata.Rad))),'k-.','LineWidth',2);
    xlabel('Radius [\mum]');
    ylabel('Occurrences [%]');
    title('Radius');
    % ISTOGRAMMA: Eccentricità
    subplot(2,2,4);
    histogram(mvn.branchdata.Eccent(:),'FaceColor',[234,179,10]./255,...
        'Normalization','probability');
    axis tight
    title('Eccentricity');
    xline(mean(mvn.branchdata.Eccent(:)),'k-.','LineWidth',2);
    xlabel('Eccentricity [adim.]');
    ylabel('Occurrences [%]');
end
%% TORTUOSITY
% Viene riplottato lo scheletro, ma questa volta il colore di ogni branch
% dipende dalla tortuosità
if tortuosity
    figure('Name','Tortuosity branch-wise')
    cm = colormap(jet);
    for b = 1:tot_branches
        currcolor_idx = 1+round((mvn.branchdata.Tort_new(b)-min(mvn.branchdata.Tort_new(:)))...
            /max(mvn.branchdata.Tort(:))*255);
        plot3(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, mvn.branchdata.zPath{b},...
            'Color', cm(currcolor_idx,:),'LineWidth',3);
        hold on
        scatter3(mvn.branchdata.From(b,1), mvn.branchdata.From(b,2), mvn.branchdata.From(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
        scatter3(mvn.branchdata.To(b,1), mvn.branchdata.To(b,2), mvn.branchdata.To(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
    end
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Tort_new(:)) max(mvn.branchdata.Tort_new(:))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    colorbar
    title('Tortuosity');
    daspect([1 1 1]);
    view(3);
end
%% RADIUS
if radius
    figure('Name','Radius branch-wise')
    cool = colormap(jet);
    for b = 1:tot_branches
        currcolor_idx = ceil(mvn.branchdata.Rad(b)/max(mvn.branchdata.Rad(:))*256);
        plot3(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, mvn.branchdata.zPath{b},...
            'Color', cool(currcolor_idx,:),'LineWidth',3);
        hold on
        scatter3(mvn.branchdata.From(b,1), mvn.branchdata.From(b,2), mvn.branchdata.From(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
        scatter3(mvn.branchdata.To(b,1), mvn.branchdata.To(b,2), mvn.branchdata.To(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
    end
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Rad(:)) max(mvn.branchdata.Rad(:))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    colorbar;
    title('Radius');
    daspect([1 1 1]);
    view(3);
end
%% LENGTH
if lengthh
    figure('Name','Length branch-wise');
    cool = colormap(jet);
    for b = 1:tot_branches
        currcolor_idx = ceil(mvn.branchdata.Len(b)/max(mvn.branchdata.Len(:))*256);
        plot3(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, mvn.branchdata.zPath{b},...
            'Color', cool(currcolor_idx,:),'LineWidth',3);
        hold on
        scatter3(mvn.branchdata.From(b,1), mvn.branchdata.From(b,2), mvn.branchdata.From(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
        scatter3(mvn.branchdata.To(b,1), mvn.branchdata.To(b,2), mvn.branchdata.To(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
    end
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Len(:)) max(mvn.branchdata.Len(:))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    colorbar;
    title('Length');
    daspect([1 1 1]);
    view(3);
end
%% ECCENTRICITY
if eccentricity
    figure('Name','Eccenticity branch-wise');
    cool = colormap(jet);
    for b = 1:tot_branches
        currcolor_idx = ceil(mvn.branchdata.Eccent(b)/max(mvn.branchdata.Eccent(:))*256);
        plot3(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, mvn.branchdata.zPath{b},...
            'Color', cool(currcolor_idx,:),'LineWidth',3);
        hold on
        scatter3(mvn.branchdata.From(b,1), mvn.branchdata.From(b,2), mvn.branchdata.From(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
        scatter3(mvn.branchdata.To(b,1), mvn.branchdata.To(b,2), mvn.branchdata.To(b,3),'o',...
            'filled','MarkerFaceColor', [1 1 1]);
    end
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Eccent(:)) max(mvn.branchdata.Eccent(:))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    colorbar;
    title('Eccentricity');
    daspect([1 1 1]);
    view(3);
end
%% PTS CLASSIFICATION
if pts_classification
    figure('Name','.pts classification');
    for i=1:tot_branches
        plot3(mvn.branchdata.xInt(i,:), mvn.branchdata.yInt(i,:), mvn.branchdata.zInt(i,:),...
            'LineWidth',2);
        hold on
        text(mvn.branchdata.xInt(i,round(end/2)), mvn.branchdata.yInt(i,round(end/2)),...
            mvn.branchdata.zInt(i,round(end/2)),string(i),'FontSize',6);
        scatter3(mvn.branchdata.From(i,1),mvn.branchdata.From(i,2),mvn.branchdata.From(i,3),'o');
        text(mvn.branchdata.From(i,1),mvn.branchdata.From(i,2),mvn.branchdata.From(i,3),...
            mvn.branchdata.CatFr(i),'FontSize',10);
        scatter3(mvn.branchdata.To(i,1),mvn.branchdata.To(i,2),mvn.branchdata.To(i,3),'o');
        text(mvn.branchdata.To(i,1),mvn.branchdata.To(i,2),mvn.branchdata.To(i,3),...
            mvn.branchdata.CatTo(i),'FontSize',10);
    end
    daspect([1 1 1]);
    title('Node classification ');
end