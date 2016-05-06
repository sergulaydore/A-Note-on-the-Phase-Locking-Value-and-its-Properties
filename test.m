
clc;
clear all;
close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\august\CircStat']));

theta = pi/4;

kappaVec = [0 1 3 15];
nVec = [5:5:100];

for idxKappa = 1:length(kappaVec)
    for idxN = 1:length(nVec)
        kappa = kappaVec(idxKappa);
        n = nVec(idxN);
        
        alpha = circ_vmrnd(theta, kappa, n);
        kappaEst = circ_kappa(alpha);
%         figure; circ_plot(alpha);

% Population PLV

popPLVTrue(idxKappa,idxN) = besseli(1,kappa)/besseli(0,kappa);
popPLVEst(idxKappa,idxN) = besseli(1,kappaEst)/besseli(0,kappaEst);
samplePLV(idxKappa,idxN) = abs(mean(exp(sqrt(-1)*alpha)));

    end
end

lw = 5; fs = 15;
figure;
plot(nVec,popPLVTrue,'--','LineWidth',lw);
% figure;
% plot(nVec,abs(popPLVEst-popPLVTrue));
% 
% figure;
 hold on
%  plot(nVec,abs(samplePLV-popPLVTrue));
plot(nVec,samplePLV,'LineWidth',lw);
xlabel('Number of Observations','Fontsize',fs)
title('Sample PLV','Fontsize',fs)
% title('$I_1(\hat{\kappa}) / I_0(\hat{\kappa})$','Interpreter',...
% 'latex','Fontsize',fs)


        
