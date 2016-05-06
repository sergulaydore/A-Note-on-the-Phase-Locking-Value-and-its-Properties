clc;
clear all;
close all;

k11 = 1; k22 = 1;

k12 = [0.001:0.001:0.999];
D = k12.^2/(k11*k22);
w = k12.^2./(2*k11*k22-k12.^2);
PLVcg_par =  pi/(sqrt(2)).*(1-1./D).* ...
    [ w.^(3/2).*hypergeom([3/4 5/4],1, w.^2) + ...
    3/4*w.^(5/2).* hypergeom([5/4 7/4],2, w.^2) ];
PLVcircgauss = abs(PLVcg_par);

figure;
lw = 3;
plot(k12,PLVcircgauss,'LineWidth',lw);
xlabel('Amplitude of cross-correlation','FontSize',15);
ylabel('Gaussian PLV','FontSize',15);

% ylabel('$\mbox{PLV}_{\mbox{circgauss}}$','Interpreter','latex','FontSize',16)
set(gca,'Ytick',[0:.2:1]); set(gca,'YTickLabel',[0:.2:1]);