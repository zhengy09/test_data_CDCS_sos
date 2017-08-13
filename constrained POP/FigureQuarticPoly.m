
clc;clear;close all;

load QuarticPoly.mat

Lwidth = 1.2;
Msize  = 7;
figure;
TimeAver = TimeAver.*100;
Index = 1:9;
h1 = loglog(N(Index),TimeAver(Index,1),'b>','markersize',Msize,'MarkerFaceColor','b'); hold on; 
     loglog(N(Index),TimeAver(Index,1),'b','linewidth',Lwidth);
h2 = loglog(N(Index),TimeAver(Index,2),'m<','markersize',Msize','MarkerFaceColor','m'); hold on; 
    loglog(N(Index),TimeAver(Index,2),'m','linewidth',Lwidth);
h3 = loglog(N(Index),TimeAver(Index,3),'ks','markersize',Msize,'MarkerFaceColor','k'); hold on; 
    loglog(N(Index),TimeAver(Index,3),'k','linewidth',Lwidth);
h = legend([h1 h2 h3],'SCS-direct','SCS-indirect','CDCS-sos','Location','NorthWest');
set(h,'box','off');
set(gcf,'Position',[250 150 450 380]);box off, grid off; %xlim([400 3.5*10^3]); ylim([1,10^4])
xlim([10,50]); ylim([0 1000]);
set(gca,'XTick',[10,15,20,30,40,50]);
xlabel('Number of variables');ylabel('CPU time per 100 iterations(s)');

load ResultLinProj.mat

Lwidth = 1.2;
Msize  = 7;
figure;
Index = 1:9;
h1 = loglog(N(Index),TimeLinear(Index,1),'b>','markersize',Msize,'MarkerFaceColor','b'); hold on; 
     loglog(N(Index),TimeLinear(Index,1),'b','linewidth',Lwidth);
h2 = loglog(N(Index),TimeLinear(Index,2),'m<','markersize',Msize','MarkerFaceColor','m'); hold on; 
    loglog(N(Index),TimeLinear(Index,2),'m','linewidth',Lwidth);
h3 = loglog(N(Index),TimeLinear(Index,3),'ks','markersize',Msize,'MarkerFaceColor','k'); hold on; 
    loglog(N(Index),TimeLinear(Index,3),'k','linewidth',Lwidth);
h = legend([h1 h2 h3],'SCS-direct','SCS-indirect','CDCS-sos','Location','NorthWest');
set(h,'box','off');
set(gcf,'Position',[250 150 450 380]);box off, grid off; %xlim([400 3.5*10^3]); ylim([1,10^4])
xlim([10,50]); ylim([0 1000]);
set(gca,'XTick',[10,15,20,30,40,50]);
xlabel('Number of variables');ylabel('CPU time per 100 iterations(s)');









