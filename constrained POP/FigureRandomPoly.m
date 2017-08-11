
clc;clear;close all;

load RandomPoly.mat

Lwidth = 1.2;
Msize  = 7;
figure;
Index = 1:6;
h1 = loglog(N(Index),TimeTotal(Index,1),'go','markersize',Msize,'MarkerFaceColor','g'); hold on; 
     loglog(N(Index),TimeTotal(Index,1),'g','linewidth',Lwidth);
Index = 1:7;
h2 = loglog(N(Index),TimeTotal(Index,2),'rp','markersize',Msize','MarkerFaceColor','r'); hold on; 
    loglog(N(Index),TimeTotal(Index,2),'r','linewidth',Lwidth);
h3 = loglog(N(Index),TimeTotal(Index,3),'cd','markersize',Msize,'MarkerFaceColor','c'); hold on; 
    loglog(N(Index),TimeTotal(Index,3),'c','linewidth',Lwidth);
h4 = loglog(N(Index),TimeTotal(Index,4),'yh','markersize',Msize,'MarkerFaceColor','y'); hold on; 
    loglog(N(Index),TimeTotal(Index,4),'y','linewidth',Lwidth);

    Index = 1:10;
h5 = loglog(N(Index),TimeTotal(Index,5),'b>','markersize',Msize,'MarkerFaceColor','b'); hold on; 
     loglog(N(Index),TimeTotal(Index,5),'b','linewidth',Lwidth);
h6 = loglog(N(Index),TimeTotal(Index,6),'m<','markersize',Msize','MarkerFaceColor','m'); hold on; 
    loglog(N(Index),TimeTotal(Index,6),'m','linewidth',Lwidth);
h7 = loglog(N(Index),TimeTotal(Index,7),'ks','markersize',Msize,'MarkerFaceColor','k'); hold on; 
    loglog(N(Index),TimeTotal(Index,7),'k','linewidth',Lwidth);

    h = legend([h1 h2 h3 h4 h5 h6 h7],'SeDuMi','SDPT3','SDPA','CSDPT',...
             'SCS-direct','SCS-indirect','CDCS-sos','location','east');
set(h,'box','off');
set(gcf,'Position',[250 150 500 400]);box off, grid off; %xlim([400 3.5*10^3]); ylim([1,10^4])
xlim([10,40]); ylim([0 1000]);
set(gca,'XTick',[10,15,20,30,40]);
xlabel('Number of variables');ylabel('Time (s)');

%print(gcf,['F:\DistributedChordalSDP\Notes\SOS-partial orthogonality\Automatica\Figures\TimeConsOpt'],'-painters','-depsc2','-r600')


Lwidth = 1.2;
Msize  = 7;
figure;
TimeAver = TimeAver.*100;
    Index = 1:10;
h5 = loglog(N(Index),TimeAver(Index,1),'b>','markersize',Msize,'MarkerFaceColor','b'); hold on; 
     loglog(N(Index),TimeAver(Index,1),'b','linewidth',Lwidth);
h6 = loglog(N(Index),TimeAver(Index,2),'m<','markersize',Msize','MarkerFaceColor','m'); hold on; 
    loglog(N(Index),TimeAver(Index,2),'m','linewidth',Lwidth);
h7 = loglog(N(Index),TimeAver(Index,3),'ks','markersize',Msize,'MarkerFaceColor','k'); hold on; 
    loglog(N(Index),TimeAver(Index,3),'k','linewidth',Lwidth);

    h = legend([h5 h6 h7],'SCS-direct','SCS-indirect','CDCS-sos','location','northwest');
set(h,'box','off');
set(gcf,'Position',[250 150 450 400]);box off, grid off; %xlim([400 3.5*10^3]); ylim([1,10^4])
xlim([10,40]); ylim([0 100]);
set(gca,'XTick',[10,15,20,30,40]);
xlabel('Number of variables');ylabel('CPU Time (s) per 100 iterations(s)');

print(gcf,['F:\DistributedChordalSDP\Notes\SOS-partial orthogonality\Automatica\Figures\TimeConsOptAver'],'-painters','-depsc2','-r600')
