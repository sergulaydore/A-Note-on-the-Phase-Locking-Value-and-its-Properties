clc;
 clear all; close all;
 
 valuesbiPLV
 
biPLV = valuesAll;
 
 valuesPLV;
PLV = valuesAll;
  
  figure;
  plot(PLV,biPLV,'*');
  axis([0 1 0 1])
  fs = 14;
%   xlabel('$\hat{\mbox{PLV}}_{\mbox{vonmises}}$','Interpreter','latex','Fontsize',fs);
%   ylabel('$\hat{\mbox{PLV}}_{\mbox{circgauss}}$','Interpreter','latex','Fontsize',fs);
xlabel('Von Mises based PLV')
xlabel('Circularly symetric Gaussian based PLV')