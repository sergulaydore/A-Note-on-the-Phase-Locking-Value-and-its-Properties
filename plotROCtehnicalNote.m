clc; clear all; close all;

load singleROCtechnicalNote.mat
lw = 3; fontsize = 15;
figure;
plot(FPcg,TPcg,'lineWidth',lw);
hold on;
plot(FPsample,TPsample,'lineWidth',lw,'Color','r');
set(gca,'Ytick',[0:0.1:1]);
set(gca,'YTickLabel',[0:0.1:1],'FontSize',fontsize);
set(gca,'Xtick',[0:0.1:1]);
set(gca,'XTickLabel',[0:0.1:1],'FontSize',fontsize);
xlabel('False positive rate','FontSize',20)
ylabel('True positive rate','FontSize',20);
h = legend('$\hat{\mbox{PLV}}_{\mbox{circgauss}}$','$\hat{\mbox{PLV}}_{\mbox{sample}}$');
 set(h,'Interpreter','latex','fontsize',13);

axis([0 .5 0.5 1])


load ROCtechnicalNote.mat

% auc : nCoupling x nSigma x nLength
% default values for
% coupling = 0.15, sigma = 0.5, n = 5000
% function of epsilon
cg = squeeze(auc_cg(:,3,3));
sample = squeeze(auc_sample(:,3,3));
aucCompare = [cg  sample];
figure; bar(aucCompare,'grouped');
h = legend('$\hat{\mbox{PLV}}_{\mbox{circgauss}}$','$\hat{\mbox{PLV}}_{\mbox{sample}}$');
 set(h,'Interpreter','latex','fontsize',13);
 axis([0 6 0.5 1.2])

set(gca,'XTickLabel',couplingVec,'FontSize',fontsize);
set(gca,'Ytick',[0:0.5:1]);
set(gca,'YTickLabel',[0:0.5:1],'FontSize',fontsize);
xlabel('\epsilon','FontSize',20)
ylabel('AUC')

% function of sigma
cg = squeeze(auc_cg(3,:,3));
sample = squeeze(auc_sample(3,:,3));
aucCompare = [cg'  sample'];
figure; bar(aucCompare,'grouped');
h = legend('$\hat{\mbox{PLV}}_{\mbox{circgauss}}$','$\hat{\mbox{PLV}}_{\mbox{sample}}$');
 set(h,'Interpreter','latex','fontsize',13);
axis([0 6 0.5 1.2])
fontsize = 15;
set(gca,'XTickLabel',sigmaVec,'FontSize',fontsize);
set(gca,'Ytick',[0:0.5:1]);
set(gca,'YTickLabel',[0:0.5:1],'FontSize',fontsize);
xlabel('\sigma','FontSize',20)
ylabel('AUC')

% function of sigma
cg = squeeze(auc_cg(3,3,:));
sample = squeeze(auc_sample(3,3,:));
aucCompare = [cg  sample];
figure; bar(aucCompare,'grouped');
h = legend('$\hat{\mbox{PLV}}_{\mbox{circgauss}}$','$\hat{\mbox{PLV}}_{\mbox{sample}}$');
 set(h,'Interpreter','latex','fontsize',13);
axis([0 6 0.5 1.2])
fontsize = 15;
set(gca,'XTickLabel',nVec,'FontSize',fontsize);
set(gca,'Ytick',[0:0.5:1]);
set(gca,'YTickLabel',[0:0.5:1],'FontSize',fontsize);
xlabel('L','FontSize',20)
ylabel('AUC')

