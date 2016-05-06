clc;
clear all;
close all;

mu = 0;
k = 1;
kappa12vec = [0 0.01 0.1  0.4 0.7 0.9];   %[0 0.5 1 2 4 8];
Qdiff = [-pi:2*pi/999:pi];
cc=hsv(length(kappa12vec));
figure
for kappa12index = 1:length(kappa12vec)
    kappa12 = -kappa12vec(kappa12index);
    a = kappa12/k *cos(Qdiff - mu);  
    q = (1./(1-a.^2)-a.*(acos(a))./sqrt((1-a.^2).^3));
    q = q/trapz(Qdiff,q);

%     [p alpha] = circ_vmpdf(Qdiff, mu, kappa);
     h = plot(Qdiff,q);
     axis([-4, 4,0,1.2]);
    set(h, 'LineStyle', '-', 'LineWidth', 2.0, 'Color',cc(kappa12index,:));
    set(gca,'Xtick',[-pi:pi:pi]); set(gca,'XTickLabel',{'p','0','p'},'FontName','Symbol');
%     set (gca,'FontName','Symbol');
% set(gca,'Ytick',[0:0.1:.4]); set(gca,'YTickLabel',[0:0.1:.4])
    hold on
    s{kappa12index}=sprintf('%s%g','|R_{12}|=',-kappa12);
end

legend(s,'FontName','Times'); xlabel('\phi','FontSize',16);
ylabel('p(\phi)','FontSize',16,'FontName','Letter')


 
    