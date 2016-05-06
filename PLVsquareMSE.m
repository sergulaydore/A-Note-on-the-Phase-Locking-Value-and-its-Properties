
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
        
        for idxMC = 1:MC
            alpha = circ_vmrnd(theta, kappa, n);
            %             kappaEst = circ_kappa(alpha);
            %         figure; circ_plot(alpha);
            popPLVTrueSq(idxKappa,idxN,idxMC) = (besseli(1,kappa)/besseli(0,kappa))^2;
            samplePLVsq(idxKappa,idxN,idxMC) = (abs(mean(exp(sqrt(-1)*alpha))))^2;
            samplePLVsqUb(idxKappa,idxN,idxMC)= (n*samplePLVsq(idxKappa,idxN,idxMC)-1)/(n-1);
            % Population PLV
        end
        
        % popPLVEstSq(idxKappa,idxN) = (besseli(1,kappaEst)/besseli(0,kappaEst))^2;
        
    end
end

msePLVSq = mean((samplePLVsq-popPLVTrueSq).^2,3);
msePLVSqUb = mean((samplePLVsqUb-popPLVTrueSq).^2,3);

lw = 5; fs = 15;
figure;
plot(nVec,msePLVSq,'LineWidth',lw);
title('$E[|\hat{PLV}^2 -(I_1({\kappa}) / I_0({\kappa}))^2|^2]$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 0 .1])

figure;
plot(nVec,msePLVSqUb,'LineWidth',lw);
title('$E[|\hat{PLV}_{ub}^2 -(I_1({\kappa}) / I_0({\kappa}))^2|^2]$',...
    'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 0 .1])
