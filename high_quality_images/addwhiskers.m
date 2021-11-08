function addwhiskers(m,s,varargin)
% ADDWHISKERS (mean, std) adds whiskers for the mean and standard
%   devitation above a plot. The whisker is displayed as
%   |------x------|
%   where the 'x' represents the mean and the '|' are in the position of
%   'mean+-std'.
%
%   Additional parametes:
%   'LegendVisibility': Allows the plot to be added to the legend
%       'on'/'off'  (default: 'on')
%   'Color': Color of the whisker
%       [rgb triplet] (default: black [0,0,0])
%   'Marker': Type of marker used for the mean
%       'o','x','*',... (default: 'x')

legendyes = 'on';
color= [0,0,0];
marker = 'x';
for i=1:numel(varargin)
    if ischar(varargin{i})
    switch varargin{i}
        case 'LegendVisibility'
            legendyes = varargin{i+1};
        case 'Color'
            color = varargin{i+1};
        case 'Marker'
            marker = varargin{i+1};
    end
    end
    
end
yl = ylim;
hg = 0.015*yl(2); h = yl(2)*1.05;
ylim([yl(1) h+3*hg]);
hold on
scatter(m,h,marker,'MarkerFaceColor',color,'MarkerEdgeColor',color,'LineWidth',1,'HandleVisibility',legendyes);
plot([m+s m+s],[h-hg h+hg],'Color',color,'HandleVisibility','off','LineWidth',1);
plot([m-s m-s],[h-hg h+hg],'Color',color,'HandleVisibility','off','LineWidth',1);
plot([m-s m+s],[h h],'Color',color,'HandleVisibility','off','LineWidth',1);
hold off
end

