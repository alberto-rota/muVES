% =============================== µVES ================================== %
% Copyrights © 2021     Alberto Rota, Luca Possenti
%
% For informations please contact:
%   alberto2.rota@mail.polimi.it
%   or alberto_rota@outlook.com
%   luca.possenti@mail.polimi.it
% ========================================================================%
%% VISUALIZATION OF 2D MICROVASCULAR NETWORK
tot_branches = size(mvn.branchdata,1);
[h,w,l] = size(mvn.bw);
[x_sk,y_sk]=ind2sub([h,w],find(mvn.skel.sk));
[x_ep,y_ep]=ind2sub([h,w],find(mvn.skel.bp));
[x_bp,y_bp]=ind2sub([h,w],find(mvn.skel.ep));
%% RAW IMAGE
if flat_img
    figure('Name','RAW Microscopy image');
    imshow(mvn.raw);
    title('RAW image from microscopy');
end
%% SEGMENTATION AND SKELETON
if segmentation || skeleton
    figure('Name','Segmentation and Skeleton');
    imshow(cat(3,mvn.bw*0.7, zeros(size(mvn.bw)),zeros(size(mvn.bw))));
    % set(gca,'Color',[0.2 0.2 0.2]);
    set(gca,'YDir','normal')
    hold on
    for i=1:size(mvn.branchdata,1)
        fnplt(mvn.branchdata.Interp{i},'y',1.5);
    end
    from = mvn.branchdata.From; to = mvn.branchdata.To;
    scatter(from(:,1),from(:,2),15,'filled','CData',[87 192 255]./255);
    scatter(to(:,1),to(:,2),15,'filled','CData',[87 192 255]./255);
    axis equal
    title('Segmentation and Skeletonization');
end
%% INTERPOLATION
if interpolation
    figure('Name','Interpolated skeleton');
    hold on;
    set(gca,'Color',[0.2 0.2 0.2]);
    title('Interpolated Skeleton');
    for i=1:tot_branches
        fnplt(mvn.branchdata.Interp{i},'y',2);
    end
    from = mvn.branchdata.From; to = mvn.branchdata.To;
    scatter(from(:,1),from(:,2),'filled','CData',[87 192 255]./255);
    scatter(to(:,1),to(:,2),'filled','CData',[255 0 0]./255);
end
%% GRAPH
if graph
    figure('Name','Graph');
    ph = plot(mvn.skel.graph,'XData',mvn.skel.graph.Nodes.x,'YData',mvn.skel.graph.Nodes.y);
    title('Graph');
    axis equal;
    for sn=1:max(mvn.branchdata.subN)
        highlight(ph,mvn.skel.graph.Nodes.subN == sn,'NodeColor',abs(rand([1 3])),...
            'MarkerSize',4);
    end
end
%% PTS CLASSIFICATION
if pts_classification
    % PLOTTAGGIO DELLA RETE CON CLASSIFICAZIONE DEI NODI
    figure('Name','.pts classification');
    for i=1:tot_branches
        plot(mvn.branchdata.xInt(i,:), mvn.branchdata.yInt(i,:),...
            'LineWidth',2);
        hold on
        text(mvn.branchdata.xInt(i,round(end/2)), mvn.branchdata.yInt(i,round(end/2)),...
            string(i),'FontSize',6);
                scatter(mvn.branchdata.From(i,1),mvn.branchdata.From(i,2),'o');
                text(mvn.branchdata.From(i,1),mvn.branchdata.From(i,2),...
                    mvn.branchdata.CatFr(i),'FontSize',10);
                scatter(mvn.branchdata.To(i,1),mvn.branchdata.To(i,2),'o');
                text(mvn.branchdata.To(i,1),mvn.branchdata.To(i,2),...
                    mvn.branchdata.CatTo(i),'FontSize',10);
    end
    title('Node classification');
end

%% HISTOGRAMS
if histograms
    % ISTOGRAMMA: Lenhezze
    figure('Name','Histogram for vessel metrics');
    subplot(1,3,1);
    histogram(mvn.branchdata.Len(:),'FaceColor','r','Normalization','probability');
    xline(mean(mvn.branchdata.Len(:)),'k-.','LineWidth',2);
    title('Length Distribution');
    xlabel('Length [\mum]');
    ylabel('Occurrences [%]');
    % ISTOGRAMMA: Tortuosità
    subplot(1,3,2);
    histogram(mvn.branchdata.Tort(:),'FaceColor','b','Normalization','probability');
    xline(mean(mvn.branchdata.Tort(:)),'k-.','LineWidth',2);
    xlabel('Tortuosity [adim.]');
    ylabel('Occurrences [%]');
    title('Tortuosity Distribution');
    % ISTOGRAMMA: Radgio
    subplot(1,3,3);
    histogram(mvn.branchdata.Rad(~isnan(mvn.branchdata.Rad)),'FaceColor','g',...
        'Normalization','probability');
        xline(median(mvn.branchdata.Rad(~isnan(mvn.branchdata.Rad))),'k-.','LineWidth',2);
    xlabel('Radius [\mum]');
    ylabel('Occurrences [%]');
    title('Radius Distribution');
end

clear num_orizz num_diagp num_diag3 distsq_to_next b n d len midp ...
    xslice yslice geo_x geo_y vesslice_x vesslice_y r_x r_y

%% TORTUOSITY
from = mvn.branchdata.From; to = mvn.branchdata.To;
if tortuosity
    figure('Name','Tortuosity')
    imshow(cat(3, mvn.bw*0.5,mvn.bw*0.5,mvn.bw*0.5));
    hold on
    cm = colormap(jet);
    for b = 1:tot_branches
        currcolor_idx = 1+round((mvn.branchdata.Tort(b)-min(mvn.branchdata.Tort(:)))...
            /max(mvn.branchdata.Tort(:))*256);
        plot(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, ...
            'Color', cm(currcolor_idx,:),'LineWidth',3);
    end
    scatter(from(:,1),from(:,2),15,'filled','CData',[1 1 1]);
    scatter(to(:,1),to(:,2),15,'filled','CData',[1 1 1]);
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Tort(:)) max(mvn.branchdata.Tort(:))]);
    set(gca,'Color',[0,0,0]);
    title('Tortuosity Branchmap');
end
%% RADIUS
if radius
    figure('Name','Radius')
    imshow(cat(3, mvn.bw*0.5,mvn.bw*0.5,mvn.bw*0.5));
    hold on
    cool = colormap(jet);
    outl = isoutlier(mvn.branchdata.Rad);
    for b = 1:tot_branches
        if any(b==find(outl))
            plot(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b},...
                'Color',[0.2 0.2 0.2],'LineWidth',3);
        else
            currcolor_idx = ceil((mvn.branchdata.Rad(b)+eps)/max(mvn.branchdata.Rad(~outl))*255);
            plot(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b},...
                'Color', cool(currcolor_idx,:),'LineWidth',3);
            hold on
        end
    end
    scatter(from(:,1),from(:,2),15,'filled','CData',[1 1 1]);
    scatter(to(:,1),to(:,2),15,'filled','CData',[1 1 1]);
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Rad(:)) max(mvn.branchdata.Rad(~outl))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    title('Radius Branchmap');
    daspect([1 1 1]);
end
%% LENGTH
if lengthh
    figure('Name','Length');
    cool = colormap(jet);
    imshow(cat(3, mvn.bw*0.5,mvn.bw*0.5,mvn.bw*0.5));
    hold on
    for b = 1:tot_branches
        currcolor_idx = ceil(mvn.branchdata.Len(b)/max(mvn.branchdata.Len(:))*256);
        plot(mvn.branchdata.xPath{b}, mvn.branchdata.yPath{b}, ...
            'Color', cool(currcolor_idx,:),'LineWidth',3);
    end
    scatter(from(:,1),from(:,2),15,'filled','CData',[1 1 1]);
    scatter(to(:,1),to(:,2),15,'filled','CData',[1 1 1]);
    hold off
    set(gca,'CLim',[min(mvn.branchdata.Len(:)) max(mvn.branchdata.Len(:))]);
    set(gca,'Color',[0.2 0.2 0.2]);
    colorbar;
    title('Length Branchmap');
end

clear b cm