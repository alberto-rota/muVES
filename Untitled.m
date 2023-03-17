close all; clear; interpreterlatex;
r = 0:25;
p = 1-exp(-0.2*r)+0.05*randn(1,numel(r)); 
p(p>=1)=1; p(p<=0)=0;
dp = [0,diff(p)];

wp = 1;
wdp = 1;

wfigure;
subplot(211); 
yyaxis left;
plot(r,p,'-o',LineWidth=1.2,MarkerSize=2.5);
grid on; hold on; 
xlabel("Repetition"); ylabel("Performance"); 
ylim([-0.5,1.5]);
a = 1-p.^(wp);
yyaxis right;  ylabel("Assistance"); 
ylim([-0.5,1.5]);
plot(r,a,'-o',LineWidth=0.5,MarkerSize=2.5);
subplot(212); 
yyaxis left;
plot(r,dp,'-o',LineWidth=1.2,MarkerSize=2.5);
a = 1-p.^(wp) - wdp*dp;
a(a<=0)=0;
grid on; hold on; 
xlabel("Repetition"); ylabel("$\delta$Performance");
ylim([-0.5,1.5]);
yyaxis right
plot(r,a,'-o',LineWidth=1.2,MarkerSize=2.5);
ylabel("Assistance");
ylim([-0.5,1.5]);