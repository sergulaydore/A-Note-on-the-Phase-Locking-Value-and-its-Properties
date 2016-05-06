
clc; clear all; close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\august\CircStat']));

theta = pi/4;

kappaVec = [0 1 3 15];
nVec = [5:5:100];
MC = 100000;

% for idxKappa = 1:length(kappaVec)
%     for idxN = 1:length(nVec)
%         kappa = kappaVec(idxKappa);
%         n = nVec(idxN);
%         popPLVTrueSq(idxKappa,idxN) = (besseli(1,kappa)/besseli(0,kappa))^2;
%         popPLVTrue(idxKappa,idxN) = (besseli(1,kappa)/besseli(0,kappa));
%         for idxMC = 1:MC
%             alpha = circ_vmrnd(theta, kappa, n);
% %             kappaEst = circ_kappa(alpha);
%             %         figure; circ_plot(alpha);
%             singlePLV(idxMC) = (abs(mean(exp(sqrt(-1)*alpha))));
%             singlePLVsq(idxMC) = (abs(mean(exp(sqrt(-1)*alpha))))^2;
%             singlePLVsqUb(idxMC)= (n*singlePLVsq(idxMC)-1)/(n-1);
%             % Population PLV
%         end
%         
%         % popPLVEstSq(idxKappa,idxN) = (besseli(1,kappaEst)/besseli(0,kappaEst))^2;
%         samplePLVmean(idxKappa,idxN) = mean(singlePLV);
%         samplePLVSqMean(idxKappa,idxN) = mean(singlePLVsq);
%         samplePLVSqUbMean(idxKappa,idxN) = mean(singlePLVsqUb);
%         
%         samplePLVvar(idxKappa,idxN) = var(singlePLV);
%         samplePLVSqVar(idxKappa,idxN) = var(singlePLVsq);
%         samplePLVSqUbVar(idxKappa,idxN) = var(singlePLVsqUb);
%     end
% end
% 
lw = 3; fs = 15;
% 
% %% Mean Plots
% 
% save samplePLVmean samplePLVmean
% save samplePLVSqMean samplePLVSqMean
% save samplePLVSqUbMean samplePLVSqUbMean
% 
% save samplePLVvar samplePLVvar
% save samplePLVSqVar samplePLVSqVar
% save samplePLVSqUbVar samplePLVSqUbVar

load samplePLVmean
load samplePLVSqMean
load samplePLVSqUbMean

load samplePLVvar
load samplePLVSqVar
load samplePLVSqUbVar

figure;
% plot(nVec,popPLVTrue,'--','LineWidth',lw); hold on;
for idx = 1:length(kappaVec)
plot(nVec,samplePLVmean,'k','LineWidth',lw);
 ylabel('$\mbox{E}\left[{PLV}_{\mbox{sample}}\right]$',...
     'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 -0.2 1.2])
end

figure;
% plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
for idx = 1: length(kappaVec)
plot(nVec,samplePLVSqUbMean(idx,:),'-b','LineWidth',lw);
% title('$E[\hat{PLV}^2_{ub\_sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
axis([0 100 -0.2 1.2])

hold on;
% plot(nVec,popPLVTrueSq,'--','LineWidth',lw); hold on;
plot(nVec,samplePLVSqMean(idx,:),'--r','LineWidth',lw);
% title('$E[\hat{PLV}^2_{sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% axis([0 100 -0.2 1.2])
h = legend('$\mbox{E} \left[\hat{\mbox{PLV}}^2_{\mbox{ub\_sample}} \right]$',...
    '$\mbox{E}\left[\hat{\mbox{PLV}}^2_{\mbox{sample}}\right]$');
set(h,'Interpreter','latex','Fontsize',13)
end


%% Variance
% figure;
% plot(nVec,log10(samplePLVvar(idx,:)),'.-k','LineWidth',lw);
% title('$Var[\hat{PLV}_{sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)

figure
for idx = 1: length(kappaVec)

%  axis([0 100 0 .09])
plot(nVec,log10(samplePLVSqUbVar(idx,:)),'-b','LineWidth',lw);
% title('$Var[\hat{PLV}^2_{ub\_sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
%  axis([0 100 0 .09])

hold on;
plot(nVec,log10(samplePLVSqVar(idx,:)),'--r','LineWidth',lw);
% title('$Var[\hat{PLV}^2_{sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
xlabel('Number of Observations','Fontsize',fs)
%  axis([0 100 0 .09])

h = legend('$\mbox{Var} \left[\hat{\mbox{PLV}}^2_{\mbox{ub\_sample}} \right]$',...
    '$\mbox{Var}\left[\hat{\mbox{PLV}}^2_{\mbox{sample}}\right]$');
set(h,'Interpreter','latex','Fontsize',fs)
end

% figure;
% plot(nVec,log10(samplePLVvar),'*-','LineWidth',lw);
% title('$Var[\hat{PLV}_{sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% %  axis([0 100 0 .09])
% 
% hold on;
% plot(nVec,log10(samplePLVSqUbVar),'--','LineWidth',lw);
% title('$Var[\hat{PLV}^2_{ub\_sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% %  axis([0 100 0 .09])
% 
% 
% plot(nVec,log10(samplePLVSqVar),'LineWidth',lw);
% title('$Var[\hat{PLV}^2_{sample]}$',...
%     'Interpreter','latex','Fontsize',fs)
% xlabel('Number of Observations','Fontsize',fs)
% %  axis([0 100 0 .09])

