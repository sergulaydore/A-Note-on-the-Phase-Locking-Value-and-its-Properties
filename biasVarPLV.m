
clc; clear all; close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\august\CircStat']));

theta = pi/4;

kappaVec = [0 1 3 15];
nVec = [5:5:100];
MC = 1000;

for idxKappa = 1:length(kappaVec)
    for idxN = 1:length(nVec)
        kappa = kappaVec(idxKappa);
         n = nVec(idxN);
%         popPLVTrueSq(idxKappa,idxN) = (besseli(1,kappa)/besseli(0,kappa))^2;
        popPLVTrue(idxKappa,idxN) = (besseli(1,kappa)/besseli(0,kappa));
        for idxMC = 1:MC
            alpha = circ_vmrnd(theta, kappa, n);
             kappaEst = circ_kappa(alpha);
            %         figure; circ_plot(alpha);
            singlePLV(idxMC) = (abs(mean(exp(sqrt(-1)*alpha))));
%             singlePLVsq(idxMC) = (abs(mean(exp(sqrt(-1)*alpha))))^2;
            singlePLVub(idxMC)= sqrt((n*singlePLVsq(idxMC)-1)/(n-1));
            singlePLVvonmises(idxMC) = (besseli(1,kappaEst)/besseli(0,kappaEst));
            % Population PLV
        end
        
        % popPLVEstSq(idxKappa,idxN) = (besseli(1,kappaEst)/besseli(0,kappaEst))^2;
        samplePLVmean(idxKappa,idxN) = mean(singlePLV);
        samplePLVubMean(idxKappa,idxN) = mean(singlePLVub);
        PLVvonmisesMean(idxKappa,idxN) = mean(singlePLVvonmises);
        
        samplePLVvar(idxKappa,idxN) = var(singlePLV);
        samplePLVubVar(idxKappa,idxN) = var(singlePLVub);
    end
end

lw = 3; fs = 15;

%% Mean Plots

figure;
% plot(nVec,popPLVTrue,'--','LineWidth',lw); hold on;
for idx = 1:length(kappaVec)
plot(nVec,samplePLVmean,'k','LineWidth',lw);
 ylabel('$\mbox{E}\left[{PLV}_{sample}\right]$',...
     'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 -0.2 1.2])
end

% figure;
% % plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
% for idx = 1: length(kappaVec)
% plot(nVec,samplePLVSqUbMean(idx,:),'-b','LineWidth',lw);
% % title('$E[\hat{PLV}^2_{ub\_sample]}$',...
% %     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% axis([0 100 -0.2 1.2])
% 
% hold on;
% % plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
% plot(nVec,samplePLVSqMean(idx,:),'--r','LineWidth',lw);
% % title('$E[\hat{PLV}^2_{sample]}$',...
% %     'Interpreter','latex','Fontsize',fs)
% % xlabel('Number of Observations','Fontsize',fs)
% % axis([0 100 -0.2 1.2])
% h = legend('$\mbox{E} \left[{PLV}^2_{ub\_sample} \right]$',...
%     '$\mbox{E}\left[{PLV}^2_{sample}\right]$');
% set(h,'Interpreter','latex')
% end
% 
% 
% %% Variance
% 
% figure
% for idx = 1: length(kappaVec)
% 
% %  axis([0 100 0 .09])
% plot(nVec,log10(samplePLVSqUbVar(idx,:)),'-b','LineWidth',lw);
% % title('$Var[\hat{PLV}^2_{ub\_sample]}$',...
% %     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% %  axis([0 100 0 .09])
% 
% hold on;
% plot(nVec,log10(samplePLVSqVar(idx,:)),'--r','LineWidth',lw);
% % title('$Var[\hat{PLV}^2_{sample]}$',...
% %     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% %  axis([0 100 0 .09])
% 
% h = legend('$\mbox{Var} \left[{PLV}^2_{ub\_sample} \right]$',...
%     '$\mbox{Var}\left[{PLV}^2_{sample}\right]$');
% set(h,'Interpreter','latex')
% end


