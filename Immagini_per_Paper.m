
%% istogramma raggi validazione 2d
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];

% !!!!!!!!!!!Cambiare questo path!!!!!!!!!!!!!!!!!
MOTHER_FOLDER = "C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice";

figure('Name','raggi_val2d','Position', [488 342 377.8000 420])
% subplot(121);
load(MOTHER_FOLDER+"\Validazione 2D\2d.mat"); 
load(MOTHER_FOLDER+"\Validazione 2D\manuale.mat"); 
load(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\dati reaver\1.mat");
% r_man = rmoutliers(r_man);
% r_man=r_man*1.59/1.98;
histogram(rmoutliers(mvn.branchdata.Rad),'FaceColor',blue,'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability','BinEdges',0:15:150);
yticklabels(yticks*100);
hold on
histogram(r_man,'FaceColor',[0.8500, 0.3250, 0.0980],'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability','BinEdges',0:15:150);
addwhiskers(mean(rmoutliers(mvn.branchdata.Rad)),std(rmoutliers(mvn.branchdata.Rad)),'Color',blue,'LegendVisibility','off');
addwhiskers(mean(r_man),std(r_man),'Color',[0.8500, 0.3250, 0.0980],'LegendVisibility','off');
%xline(mvn.mRad,'Color',blue, 'LineWidth',1.5);
hold on
scatter(metrics.meanVesselDiam/2*1.59,0.62,[],green,'x', 'LineWidth',1.5);
% [0, 0.4470, 0.7410]
%xline(mean(r_man),'Color',orange, 'LineWidth',1.5);
ylim([0 0.65]);
grid on
xlabel('Radius [\mum]')
ylabel('Occurrences [%]')
legend('Mean Radius with \muVES','Mean manual Radius','Mean Radius with REAVER','Location','northoutside');
axis square
%% istogramma cumulative length
figure('Name','length_val2d','Position', [488 342 377.8000 420])
[hl,ed] = histcounts(mvn.branchdata.Len,100);
hl = cumsum(hl.*ed(2:end));
sum(mvn.branchdata.Len);
stairs(ed(2:end),hl,'LineWidth',1,'Color',[0, 0.4470, 0.7410]); 
hold on
[hl,ed] = histcounts(l_man,100);
hl = cumsum(hl.*ed(2:end));
sum(mvn.branchdata.Len*1.7);
stairs(ed(2:end),hl,'LineWidth',1,'Color',orange);
grid on
ylabel('Cumulative Length [mm]');
xlabel('Branch Length [\mum]');
yticklabels(yticks/1000);

yline(metrics.vesselLength_CORR*1.59,'Color',green, 'LineWidth',1.5);
%  xline(sum(l_man),'Color',orange, 'LineWidth',1.5);
% xline(metrics.vesselLength*1.59,'Color',green,'LineWidth',1.5);
legend('\muVES','Manual Analysis','REAVER total length','Location','northoutside');
axis square
%% rete2d+skeletro
load(MOTHER_FOLDER+"\Validazione 2D\2d.mat"); 
figure;
imshow(cat(3,mvn.bw*0.7, zeros(size(mvn.bw)),zeros(size(mvn.bw))));
set(gca,'Color',[0.2 0.2 0.2]);
hold on
for i=1:size(mvn.branchdata,1)
fnplt(mvn.branchdata.Interp{i},'y',2);
end
from = mvn.branchdata.From; to = mvn.branchdata.To;
scatter3(from(:,1),from(:,2),from(:,3),'filled','CData',[87 192 255]./255);
scatter3(to(:,1),to(:,2),to(:,3),'filled','CData',[255 0 0]./255);

%% rete3d_dws2+skeletro
load(MOTHER_FOLDER+"\Da 2D a 3D\3ddws.mat");
figure;
p = patch(isosurface(mvn.bw));
reducepatch(p,0.1)
daspect([1 1 1]);
view(3);
p.FaceColor = 'red'; p.EdgeColor = 'red';
p.FaceAlpha = 0.01;
p.EdgeAlpha = 0.1;
set(gca,'Color',[0.2 0.2 0.2]);
hold on
for i=1:size(mvn.branchdata,1)
fnplt(mvn.branchdata.Interp{i},'y',2);
end
from = mvn.branchdata.From; to = mvn.branchdata.To;
scatter3(from(:,1),from(:,2),from(:,3),'filled','CData',[87 192 255]./255);
scatter3(to(:,1),to(:,2),to(:,3),'filled','CData',[255 0 0]./255);
ca = gca;
ca.Box = 'on';
ca.BoxStyle = 'full';
ca.XColor = [1 1 1];
ca.YColor = [1 1 1];
ca.ZColor = [1 1 1];
ca.LineWidth = 1;
title('Downsampled');
view([-27 47]);
%% rete3d+skeletro
load(MOTHER_FOLDER+"\Da 2D a 3D\3d.mat");
figure;
p = patch(isosurface(mvn.bw));
reducepatch(p,0.1)
daspect([1 1 1]);
view(3);
p.FaceColor = 'red'; p.EdgeColor = 'red';
p.FaceAlpha = 0.01;
p.EdgeAlpha = 0.1;
set(gca,'Color',[0.2 0.2 0.2]);
hold on
for i=1:size(mvn.branchdata,1)
fnplt(mvn.branchdata.Interp{i},'y',2);
end
from = mvn.branchdata.From; to = mvn.branchdata.To;
scatter3(from(:,1),from(:,2),from(:,3),'filled','CData',[87 192 255]./255);
scatter3(to(:,1),to(:,2),to(:,3),'filled','CData',[255 0 0]./255);
ca = gca;
ca.Box = 'on';
ca.BoxStyle = 'full';
ca.XColor = [1 1 1];
ca.YColor = [1 1 1];
ca.ZColor = [1 1 1];
ca.LineWidth = 1;
title('Full');
view([-27 47]);
%% istogramma 2d - 3d
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];

figure('Name','istogrammi_2d3d3ddws','Position', [634.6000 115.4000 393.6000 372])
ax1=subplot(321); 
load(MOTHER_FOLDER+"\Da 2D a 3D\2d.mat");
histogram(rmoutliers(mvn.branchdata.Rad),'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,180,14));
addwhiskers(mean(rmoutliers(mvn.branchdata.Rad)), std(rmoutliers(mvn.branchdata.Rad)),'Color',orange);
xline(mean(rmoutliers(mvn.branchdata.Rad)),'--','Color',orange,'LineWidth',1.2);
xlabel('Radius [\mum]')
ylabel('Occurrences [%]')
% title(['R_{mean} = ' num2str(mvn.mRad) '\mum']);    
grid on
ax3=subplot(323); 
load(MOTHER_FOLDER+"\Da 2D a 3D\3ddws.mat");
histogram(mvn.branchdata.Rad,'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,180,14));
addwhiskers(mean(mvn.branchdata.Rad), std(mvn.branchdata.Rad),'Color',orange);
xline(mvn.mRad,'--','Color',orange,'LineWidth',1.2);
xlabel('Radius [\mum]')
ylabel('Occurrences [%]')
% title(['R_{mean} = ' num2str(mvn.mRad) '\mum']);
grid on
ax5=subplot(325);
load(MOTHER_FOLDER+"\Da 2D a 3D\3d.mat");
histogram(mvn.branchdata.Rad,'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,180,14));
addwhiskers(mean(mvn.branchdata.Rad), std(mvn.branchdata.Rad),'Color',orange);
xline(mvn.mRad,'--','Color',orange,'LineWidth',1.2);
xlabel('Radius [\mum]')
ylabel('Occurrences [%]')
% title(['R_{mean} = ' num2str(mvn.mRad) '\mum']);
grid on
xlim([-10 100]);
% figure
ax2=subplot(322); 
load(MOTHER_FOLDER+"\Da 2D a 3D\2d.mat");
histogram(mvn.branchdata.Len,'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,500,14));
addwhiskers(mean(rmoutliers(mvn.branchdata.Len)), std(rmoutliers(mvn.branchdata.Len)),'Color',blue);
xline(mean(rmoutliers(mvn.branchdata.Len)),'--','Color',blue,'LineWidth',1.2);
xlabel('Length [\mum]')
ylabel('Occurrences [%]')
% title(['L_{mean} = ' num2str(mvn.mLen) '\mum']);
grid on
ax4=subplot(324); 
load(MOTHER_FOLDER+"\Da 2D a 3D\3ddws.mat");
histogram(mvn.branchdata.Len,'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,500,14));
addwhiskers(mean(mvn.branchdata.Len), std(mvn.branchdata.Len),'Color',blue);
xline(mvn.mLen,'--','Color',blue,'LineWidth',1.2);
xlabel('Length [\mum]')
ylabel('Occurrences [%]')
% title(['L_{mean} = ' num2str(mvn.mLen) '\mum']);
grid on
ax6=subplot(326);
load(MOTHER_FOLDER+"\Da 2D a 3D\3d.mat");
histogram(mvn.branchdata.Len,'FaceColor',[0.7 0.7 0.7],'LineWidth',1.5,'Normalization','Probability',...
    'BinEdges',linspace(0,500,14));
addwhiskers(mean(mvn.branchdata.Len), std(mvn.branchdata.Len),'Color',blue);
xline(mvn.mLen,'--','Color',blue,'LineWidth',1.2);
xlabel('Length [\mum]')
ylabel('Occurrences [%]')
% title(['L_{mean} = ' num2str(mvn.mLen) '\mum']);
grid on
linkaxes([ax1,ax3,ax5],'x')
linkaxes([ax2,ax4,ax6],'x')
load(MOTHER_FOLDER+"\Da 2D a 3D\2d.mat");
tab(1,1:4) = [mvn.mRad mvn.mLen mvn.mTort NaN];
load(MOTHER_FOLDER+"\Da 2D a 3D\3ddws.mat");
tab(2,1:4) = [mvn.mRad mvn.mLen mvn.mTort mvn.mEcc];
load(MOTHER_FOLDER+"\Da 2D a 3D\3d.mat");
tab(3,1:4) = [mvn.mRad mvn.mLen mvn.mTort mvn.mEcc];
%% istogramma 3d DL vs AC (ramo per ramo)
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];

load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_nodws.mat");
figure('Name','istogrammi_3ddlac','Position', [634.6000 115.4000 393.6000 372])
histogram((mvn.branchdata.Len),'FaceColor',blue,'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability',...
    'BinEdges',0:500/12:500);
yticklabels(yticks*100);
hold on
% xline(mvn.mLen,'Color',blue, 'LineWidth',1.5);
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\AC_nodws.mat");
mvnac=mvn;
histogram((mvn.branchdata.Len),'FaceColor',orange,'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability',...
    'BinEdges',0:500/12:500);
addwhiskers(mean(rmoutliers(mvn.branchdata.Len)),std(mvn.branchdata.Len),'Color',orange);
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_nodws.mat");
mvndl=mvn;
addwhiskers(mean((mvn.branchdata.Len)),std(mvn.branchdata.Len),'Color',blue);
hold on
% xline(mvn.mLen,'Color',orange, 'LineWidth',1.5);
grid minor; axis square;    
legend('DL','AC');
xlabel('Length [\mum]');
ylabel('Occurrences [%]');


load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_nodws.mat");
figure('Name','istogrammi_3ddlac','Position', [634.6000 115.4000 393.6000 372])
histogram(rmoutliers(mvn.branchdata.Rad),'FaceColor',blue,'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability',...
    'BinEdges',0:5:70);
hold on
yticklabels(yticks*100);
% xline(mvn.mRad,'Color',blue, 'LineWidth',1.5);
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\AC_nodws.mat");
histogram(rmoutliers(mvn.branchdata.Rad),'FaceColor',orange,'LineWidth',1.5,'HandleVisibility','on','Normalization','Probability',...
    'BinEdges',0:5:70);
addwhiskers(mean(rmoutliers(mvn.branchdata.Rad)),std(mvn.branchdata.Rad),'Color',orange);
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_nodws.mat");
addwhiskers(mean(rmoutliers(mvn.branchdata.Rad)),std(mvn.branchdata.Rad),'Color',blue);
hold on
% xline(mvn.mRad,'Color',orange, 'LineWidth',1.5);
grid minor; axis square;
legend('DL','AC');
xlabel('Radius [\mum]');
ylabel('Occurrences [%]');
%% barplot tempi computazionali dl/ac
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\AC_nodws.mat")
mvnac= mvn;
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_nodws.mat")
mvndl= mvn;
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\AC_dws.mat")
mvnac2= mvn;
load("C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice\DL e AC Senza downsampling\DL_dws.mat")
mvndl2= mvn;
%%
mvnac = ac;
mvndl = dl;
figure('Name','barplot_chrono','Position', [169 440.2000 1.0608e+03 209.6000])
cats = categorical(["DL","AC","DL (DF=2)","AC (DF=2)"]);

times =[mvndl.info.chrono.Duration(4),...
        mvndl.info.chrono.Duration(10),...
        sum(mvndl.info.chrono.Duration([1:3,5:9,11]));
        mvnac.info.chrono.Duration(5),...
        mvnac.info.chrono.Duration(11),...
        sum(mvnac.info.chrono.Duration([1:4,6:10,12]));
        mvndl2.info.chrono.Duration(4),...
        mvndl2.info.chrono.Duration(10),...
        sum(mvndl2.info.chrono.Duration([1:3,5:9,11]));
        mvnac2.info.chrono.Duration(5),...
        mvnac2.info.chrono.Duration(11),...
        sum(mvnac2.info.chrono.Duration([1:4,6:10,12]))];

b = barh(cats,times,'stacked');
barwidth=0.5;
b(1).BarWidth = barwidth;
green = [0.4660, 0.6740, 0.1880];
b(3).FaceColor = green;
xticks(0:50:1400)
xlabel('Computational time [seconds]')
legend('Segmentaion','Morphological Measurements','Others')
grid on
%% scatterplot diagonale lunghezza e raggio(confronto valori medi) 2D [SLOW]
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
yellow = "#ffe012";
purple = "#6f3885";
r_re = zeros(1,30);
r_mu = zeros(1,30);
l_re = zeros(1,30);
l_mu = zeros(1,30);
for i=1:30
    load(strrep(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\dati reaver\?.mat","?",string(i)));
    l_re(i) = metrics.vesselLength_CORR*1.59;
    load(strrep(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\?_rp1.mat","?",string(i)));
    r_mu_RE(i) = mvn.mRad_REAVER; l_mu(i) = sum(mvn.branchdata.Len); r_re(i) = metrics.meanVesselDiam/2*1.59; 
    r_mu1(i) = mvn.mRad;
    load(strrep(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\?_rp3.mat","?",string(i)));
    r_mu3(i) = mvn.mRad;
%     load(strrep(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\?_rp5.mat","?",string(i)));
%     r_mu5(i) = mvn.mRad;
%     load(strrep(MOTHER_FOLDER+"\Confronto REAVER_Noi 30 img\?_rp10.mat","?",string(i)));
%     r_mu10(i) = mvn.mRad;
end

figure('Name','scatter_rag','Position', [488 342 377.8000 420])
plot([0,100],[0,100],'--','LineWidth',2,'Color',[0.7 0.7 0.7],'HandleVisibility','off');
hold on
scatter((r_re),(r_mu1),'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.8);
scatter((r_re),(r_mu3),'filled','MarkerFaceColor',green,'MarkerFaceAlpha',0.8);
% scatter((r_re),(r_mu5),'filled','MarkerFaceColor',yellow,'MarkerFaceAlpha',0.8);
% scatter((r_re),(r_mu10),'filled','MarkerFaceColor',purple,'MarkerFaceAlpha',0.8);
scatter((r_re),(r_mu_RE),'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.8);
% p = polyfit(sort(r_re),sort(r_mu),1);
% fplot(@(x) p(1).*x + p(2),':','Color',blue,'LineWidth',1);

% scatter(sort(r_re),sort(r_mu_RE),'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.8);
xlabel('Radius [\mum] with REAVER');
ylabel('Radius [\mum] with \muVES');
grid minor; axis square; %axis square
% xlim([min(r_re) max(r_mu)]);
% xlim([min(r_re) max(r_mu)]);
xlim([16 60]);
ylim([16 60]);
title(strcat("RMSE = ",string(mean(sqrt((r_mu1-r_re).^2))), "\mum"));
legend('Branch average [RP = 1]','Branch average [RP = 3]','Length-Weighted average',...
    'Location','southeast');

figure('Name','scatter_len','Position', [488 342 377.8000 420])
plot([0,50000],[0,50000],'--','LineWidth',2,'Color',[0.7 0.7 0.7]);
hold on
scatter((l_re/1000),(l_mu/1000),'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.8);
% p = polyfit(sort(l_re),sort(l_mu),1);
% fplot(@(x) p(1).*x + p(2),':','Color',blue,'LineWidth',1);
% xticklabels(xticks/1000);
% yticklabels(yticks/1000);
xlabel('Length [mm] with REAVER');
ylabel('Length [mm] with \muVES');
grid minor; axis square; %axis square
% xlim([min(l_re) max(l_mu)]);
% xlim([min(l_re) max(l_mu)]);
xlim([5 20]);
title(strcat("RMSE = ",string(mean(sqrt((l_mu-l_re).^2))/1000), "mm"));
%% scatterplot diagonale lunghezza e raggio(confronto valori medi) 3D [SLOW]
MOTHER_FOLDER = "C:\Users\alber\OneDrive - Politecnico di Milano\mVN_codice";
figure()
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
r_dl = zeros(1,30);
r_ac = zeros(1,30);
l_dl = zeros(1,30);
l_ac = zeros(1,30);
for i=1:30
    load(strrep(MOTHER_FOLDER+"\Confronto DL_AC 30 img\?_DL.mat","?",string(i)));
    l_dl(i) = mvn.mLen; r_dl(i) = mvn.mRad;
    load(strrep(MOTHER_FOLDER+"\Confronto DL_AC 30 img\?_AC.mat","?",string(i)));
    l_ac(i) = mvn.mLen; r_ac(i) = mvn.mRad;
end
%%
figure('Name','rag_dlac','Position', [634.6000 115.4000 393.6000 372])
plot([0,100],[0,100],'--','LineWidth',2,'Color',[0.7 0.7 0.7],'HandleVisibility','off');
hold on
scatter((r_ac),(r_dl),'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.8);
% scatter((r_re),(r_mu_RE),'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.8);
% p = polyfit(sort(r_re),sort(r_mu),1);
% fplot(@(x) p(1).*x + p(2),':','Color',blue,'LineWidth',1);

% scatter(sort(r_re),sort(r_mu_RE),'filled','MarkerFaceColor',orange,'MarkerFaceAlpha',0.8);
xlabel('Radius [\mum] after AC');
ylabel('Radius [\mum] after DL');
grid minor; axis square; %axis square
% xlim([min(r_re) max(r_mu)]);
% xlim([min(r_re) max(r_mu)]);
% xlim([16 60]);
% ylim([16 60]);
xlim([18 42]);
ylim([18 42]);
title(strcat("RMSE = ",string(mean(sqrt((r_dl-r_ac).^2))), "\mum"));
% legend('\muVES method','REAVER method');

figure('Name','len_dlac','Position', [634.6000 115.4000 393.6000 372])
plot([0,50000],[0,50000],'--','LineWidth',2,'Color',[0.7 0.7 0.7]);
hold on
scatter((l_dl),(l_ac),'filled','MarkerFaceColor',blue,'MarkerFaceAlpha',0.8);
% p = polyfit(sort(l_re),sort(l_mu),1);
% fplot(@(x) p(1).*x + p(2),':','Color',blue,'LineWidth',1);
% xticklabels(xticks/1000);
% yticklabels(yticks/1000);
xlabel('Length [\mum] after AC');
ylabel('Length [\mum] after DL');
grid minor; axis square; %axis square
% xlim([min(l_re) max(l_mu)]);
% xlim([min(l_re) max(l_mu)]);
% xlim([16 60]);
% ylim([16 60]);
xlim([50 135]);
ylim([50 135]);
title(strcat("RMSE = ",string(mean(sqrt((l_dl-l_ac).^2))), "\mum"));
%% tabella parametri dl-ac
load('3d_dl.mat');
tab(1,1:4) = [mvn.mRad mvn.mLen mvn.mTort mvn.mEcc];
load('3d_ac.mat');
tab(2,1:4) = [mvn.mRad mvn.mLen mvn.mTort mvn.mEcc];

%% confronto deeplearning
load('3d_ac.mat');
bwac = mvn.bw;
load('3d_dl.mat');
bwdl = mvn.bw;
v = zeros(size(bwac));
v(bwac & bwdl & bwac==1) = 1; %TP
v(bwac & bwdl & bwac==0) = 0; %TN
v(xor(bwac,bwdl) & bwac==0) = 3; %FP
v(xor(bwac,bwdl) & bwac==1) = 4; %FN
TP = nnz(v==1);
TN = nnz(v==0);
FP = nnz(v==3);
FN = nnz(v==4);
conf = [TP FN; FP TN];

%% tempi computazionali
load('3d_dl.mat');
tabdl = mvn.info.chrono;
load('3d_ac.mat');
tabac = mvn.info.chrono;
writetable(tabdl,'tempi.xlsx','WriteMode','append');
writetable(tabac,'tempi.xlsx','WriteMode','append');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAFICI CON CELLULA TUMORALE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANDAMENTO VALORI MEDI 
rad = [];
len = [];
for i=1:6
    path = strrep(MOTHER_FOLDER+"\Stack_tumorale_ok\Image0?_02.mat","?",string(i));
    load(path);
    rad = cat(1,rad,mvn.mRad);
    len = cat(1,len,mvn.mLen);
end
dist_to_tumoral = [0:mvn.info.pxdensity(1)*size(mvn.bw,1):mvn.info.pxdensity(1)*size(mvn.bw,1)*5]/1000;
yyaxis left
plot(dist_to_tumoral,rad);
ylabel("Mean Radius [\mum]")
yyaxis right
plot(dist_to_tumoral,len);
ylabel("Mean Length [\mum]")
xlabel("Distance to tumoral cell [mm]");


%% MONTAGGIO 6 STACKS QUADRATI IN UN UNICA IMMAGINE LUNGA
m = [];
for i=1:6
    path = strrep(MOTHER_FOLDER+"\Stack_tumorale_ok\Image0?_02.mat","?",string(i));
    load(path);
    m = cat(1,m,mat2gray(mvn.flat));
end
figure
imshow(m,[]);
imwrite(m,"C:\Users\alber\Desktop\mont.bmp");
%% CALCOLO DELLA SUPERFICIE INTERPOLANTE (fare una volta sola, richiede tempo...)
xv=[]; yv=[]; rv=[]; lv=[]; tv = []; ev = []; gv = [];
for i=1:6
    path = strrep(MOTHER_FOLDER+"\Stack_tumorale_ok\Image0?_02.mat","?",string(i));
    load(path);
    xv = cat(1,xv,(mvn.branchdata.From(:,1)+mvn.branchdata.To(:,1))/2);
    yv = cat(1,yv,(mvn.branchdata.From(:,2)+mvn.branchdata.To(:,2))/2+640*(i-1));
    rv = cat(1,rv,mvn.branchdata.Rad);
    lv = cat(1,lv,mvn.branchdata.Len);
    tv = cat(1,tv,mvn.branchdata.Tort);
    ev = cat(1,ev,mvn.branchdata.Eccent);
    gv = cat(1,gv,double(mvn.branchdata.isGood));
end
% [rv,remr] = rmoutliers(rv);
% xv = xv(not(remr)); yv = yv(not(remr)); lv=lv(not(remr));
% [lv,reml] = rmoutliers(lv);
% xv = xv(not(reml)); yv = yv(not(reml)); rv=rv(not(reml));

%% VISUALIZZAZIONE DELA SUPERFICIE INTERPOLANTE
Rint = fit([xv,yv],rv,'poly22');
Lint = fit([xv,yv],lv,'poly22');
Tint = fit([xv,yv],tv,'poly22');
Eint = fit([xv(not(isnan(ev))),yv(not(isnan(ev)))],ev(not(isnan(ev))),'poly22');
Gint = fit([xv,yv],gv,'poly22');

figure; %subplot(121); 
scatter3(xv,yv,rv,'.'); hold on; s = plot(Rint); colorbar 
s.EdgeAlpha = 0; shading interp; daspect([1,1,0.05]);  view([-62.16 22.13]) 
scatter3(640/2,640/2,0,70,'r','filled');
xlim([0,640]); ylim([0,3840]); zlim([0,60])
zlabel("Branch Radius [\mum]")
figure; %subplot(122);
scatter3(xv,yv,lv,'.'); hold on; s = plot(Lint); colorbar
s.EdgeAlpha = 0; shading interp; daspect([1,1,0.3]); view([-62.16 22.13]) 
scatter3(640/2,640/2,0,70,'r','filled');
zlabel("Branch Length [\mum]")
xlim([0,640]); ylim([0,3840]); zlim([0,250])
figure; %subplot(122);
scatter3(xv,yv,ev,'.'); hold on; s = plot(Eint);colorbar 
s.EdgeAlpha = 0; shading interp; daspect([1,1,0.003]);  view([-62.16 22.13]) 
scatter3(640/2,640/2,0,70,'r','filled');
zlabel("Branch Eccentricity [-]")
xlim([0,640]); ylim([0,3840]); %zlim([0,250])
figure; %subplot(122);
scatter3(xv,yv,tv,'.'); hold on; s = plot(Tint); colorbar
s.EdgeAlpha = 0; shading interp; daspect([1,1,0.0015]); view([-62.16 22.13]) 
scatter3(640/2,640/2,1,70,'r','filled');
zlabel("Branch Tortuosity [-]")
xlim([0,640]); ylim([0,3840]); zlim([1,2])
figure; %subplot(122);
scatter3(xv,yv,gv,'.'); hold on; s = plot(Gint); colorbar
s.EdgeAlpha = 0; shading interp; daspect([1,1,0.004]); view([-62.16 22.13]) 
scatter3(640/2,640/2,0,70,'r','filled')
zlabel("Branch 'isGood' [bool]")
xlim([0,640]); ylim([0,3840]); %zlim([0,250])