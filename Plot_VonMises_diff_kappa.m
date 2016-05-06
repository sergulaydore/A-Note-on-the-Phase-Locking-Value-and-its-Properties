clc;
clear all;
close all;
addpath(['F:\29.02.2012backup\old research\RESEARCH SPG 10\PROJECTS-June 2010'...
    '\Von mises\CircStat'])
mu = 0;

kappavec = [0 0.5 1 2 4 8];
Qdiff = [-pi:2*pi/999:pi];
cc=hsv(length(kappavec));
legendtex = [];
figure
for kappaindex = 1:length(kappavec)
    kappa = kappavec(kappaindex);
    [p alpha] = circ_vmpdf(Qdiff, mu, kappa);
    h = plot(alpha,p);
    axis([-4, 4,0,1.2]);
    set(h, 'LineStyle', '-', 'LineWidth', 2.0, 'Color',cc(kappaindex,:));
    set(gca,'Xtick',[-pi:pi:pi]); set(gca,'XTickLabel',{'p','0','p'},'FontName','Symbol');
%     set (gca,'FontName','Symbol');
% set(gca,'Ytick',[0:0.1:.4]); set(gca,'YTickLabel',[0:0.1:.4])
    hold on
    s{kappaindex}=sprintf('%s%g','\kappa=',kappa);
end


legend(s,'FontName','Times');
xlabel('\phi','FontSize',16)
ylabel('p(\phi)','FontSize',16,'FontName','Letter')


 
    