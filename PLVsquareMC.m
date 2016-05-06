
clc;
clear all;
close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\august\CircStat']));

theta = pi/4;

kappaVec = [0 1 3 15];
nVec = [5:5:100];
MC = 1000;

for idxKappa = 1:length(kappaVec)
    for idxN = 1:length(nVec)
        kappa = kappaVec(idxKappa);
        n = nVec(idxN);
        popPLVTrueSq(idxKappa,idxN) = (besseli(1,kappa)/besseli(0,kappa))^2;
        
        for idxMC = 1:MC
            alpha = circ_vmrnd(theta, kappa, n);
%             kappaEst = circ_kappa(alpha);
            %         figure; circ_plot(alpha);
            singlePLVsq(idxMC) = (abs(mean(exp(sqrt(-1)*alpha))))^2;
            singlePLVsqUb(idxMC)= (n*singlePLVsq(idxMC)-1)/(n-1);
            % Population PLV
        end
        
        % popPLVEstSq(idxKappa,idxN) = (besseli(1,kappaEst)/besseli(0,kappaEst))^2;
        samplePLVSq(idxKappa,idxN) = mean(singlePLVsq);
        samplePLVSqUb(idxKappa,idxN) = mean(singlePLVsqUb);
    end
end

lw = 5; fs = 15;
figure;
% plot(nVec,popPLVTrueSq,'--','LineWidth',lw);
% figure;
% plot(nVec,abs(popPLVEst-popPLVTrue));
%
% figure;
%  hold on
plot(nVec,abs(samplePLVSq-popPLVTrueSq),'LineWidth',lw);
% plot(nVec,samplePLVSq,'LineWidth',lw);
title('$|E[\hat{PLV}^2] -I_1({\kappa}) / I_0({\kappa})^2|$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 0 .2])
% title('Sample PLV^2','Fontsize',fs)
% title('$I_1(\hat{\kappa}) / I_0(\hat{\kappa})$','Interpreter',...
% 'latex','Fontsize',fs)

figure;
plot(nVec,abs(samplePLVSqUb-popPLVTrueSq),'LineWidth',lw);
title('$|E[\hat{PLV}^2_{ub}] -I_1({\kappa}) / I_0({\kappa})^2|$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 0 .2])

figure;
plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
plot(nVec,samplePLVSq,'LineWidth',lw);
title('$E[\hat{PLV}^2]$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 -0.2 1.2])

figure;
plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
plot(nVec,samplePLVSqUb,'LineWidth',lw);
title('$E[\hat{PLV}^2_{ub]}$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 -0.2 1.2])