clc;
clear all;
close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\august\CircStat']));

nParam1 = 1000;
nParam2 = 1000;
muVec = [-pi:2*pi/(nParam1-1):pi];
kappaVec = [0:10/(nParam2-1):10];

n = 500;
for idx1 = 1:nParam1;
    mu = muVec(idx1);
    for idx2 = 1:nParam2;

kappa =kappaVec(idx2);
alpha = circ_vmrnd(mu, kappa, n);

PLV(idx1,idx2) = abs(mean(exp(j*alpha)));
PLI(idx1,idx2) = abs(mean(sign(alpha)));
    end
end
fs=15;
figure;
imagesc(kappaVec,muVec,PLV); caxis([0 1]); colorbar;
xlabel('Concentration Parameter','Fontsize',fs);
 ylabel('Mean','Fontsize',fs)


figure;
imagesc(kappaVec,muVec,PLI); caxis([0 1]); colorbar
xlabel('Concentration Parameter','Fontsize',fs);
 ylabel('Mean','Fontsize',fs)
