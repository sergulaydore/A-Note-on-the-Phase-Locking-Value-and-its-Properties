clc;
clear all;
close all;

if ismac
    addpath('/Volumes/VENCEREMOS/29.02.2012backup/Research Spring 2012/april/paper_draft/matlab_codes/functions');
else
    saveDir = 'C:\Users\Sergul\Documents\My Dropbox\matlab_codes\roessler\3DRoesslerROC';
    addpath('C:\Users\Sergul\Documents\My Dropbox\matlab_codes\functions')
end


fontsize = 18; markersize = 14;

n = 10000;
MC = 1000;
delta_thresh = 0.005;
couplingVec = [0.01:0.01:0.3];

for idxCoupling = 2:1:length(couplingVec)
    idxCoupling
    coupling = couplingVec(idxCoupling);
[TPcg FPcg TPnp ...
    FPnp TPsample FPsample...
    TPvm FPvm...
    auc_cg(idxCoupling) auc_np(idxCoupling) auc_sample(idxCoupling) auc_vm(idxCoupling)]...
    = f_ROC(coupling, MC, n, delta_thresh);
filename = [saveDir '\ROC' num2str(idxCoupling)];
save(filename,'TPcg','FPcg','TPnp','FPnp','TPsample','FPsample',...
    'TPvm','FPvm')
end

save([saveDir '\parameters' ] ,'auc_cg','auc_np','auc_sample','auc_vm',...
    'couplingVec','MC','n','delta_thresh')


% figure; plot(FPcg0, TPcg0, '.','MarkerSize',markersize); hold on;
% plot(FPnp0, TPnp0, '.r','MarkerSize',markersize);
% plot(FPsample0, TPsample0, '.k','MarkerSize',markersize);
% plot(FPvm0, TPvm0, '.g','MarkerSize',markersize);
% xlabel('False Positive', 'FontSize', fontsize); ylabel('True Positive', 'FontSize', fontsize);
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV','Von Mises');
% title(sprintf('when coupling %.2f', coupling0))
% set(gca, 'linewidth',2, 'FontSize', fontsize);

%  %%
%
%  coupling = 0.1;
%
% % [TPcg FPcg TPnp FPnp TPsample FPsample auc_cg auc_np auc_sample] = f_ROC(coupling, MC, n, delta_thresh);
%
%
% figure; plot(FPcg, TPcg, '.','MarkerSize',markersize); hold on;
% plot(FPnp, TPnp, '.r','MarkerSize',markersize); plot(FPsample, TPsample, '.k','MarkerSize',markersize);
% xlabel('False Positive', 'FontSize', fontsize); ylabel('True Positive', 'FontSize', fontsize);
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV');
%  title(sprintf('when coupling %.2f', coupling))
% set(gca, 'linewidth',2, 'FontSize', fontsize);
%
% % linespec = 'b';
% % figure; bar([auc_cg ; auc_np; auc_sample],linespec);
% % text(1-.1,auc_cg+.03,sprintf('%.2f',auc_cg))
% % text(2-.1,auc_np+.03,sprintf('%.2f',auc_np))
% % text(3-.1,auc_sample+.03,sprintf('%.2f',auc_sample))
% % axis([0.5 3.5 0 1.2])
% %
% % figure;
% % bar([auc_cg,0,0],'b');
% % hold on
% % bar([0,auc_np,0],'r');
% % bar([0,0,auc_sample],'k')
% %%
% coupling2 = 0.2;
% % [TPcg2 FPcg2 TPnp2 FPnp2 TPsample2 FPsample2 auc_cg2 auc_np2 auc_sample2] = f_ROC(coupling, MC, n, delta_thresh);
%
% figure; plot(FPcg2, TPcg2, '.','MarkerSize',markersize); hold on;
% plot(FPnp2, TPnp2, '.r','MarkerSize',markersize);
% plot(FPsample2, TPsample2, '.k','MarkerSize',markersize);
% xlabel('False Positive Rate', 'FontSize', fontsize);
% ylabel('True Positive Rate', 'FontSize', fontsize);
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV');
% %  title(sprintf('when coupling %.2f', coupling2))
% set(gca, 'linewidth',2, 'FontSize', fontsize);
%
%
%
% %%
%
% n = 10000;
% MC = 5000;
%  delta_thresh = 0.005;
%  coupling3 = 0.3;
% % [TPcg3 FPcg3 TPnp3 FPnp3 TPsample3 FPsample3 auc_cg3 auc_np3 auc_sample3] = f_ROC(coupling, MC, n, delta_thresh);
%
%
% figure; plot(FPcg3, TPcg3, '.','MarkerSize',markersize); hold on;
% plot(FPnp3, TPnp3, '.r','MarkerSize',markersize); plot(FPsample3, TPsample3, '.k','MarkerSize',markersize);
% xlabel('False Positive', 'FontSize', fontsize); ylabel('True Positive', 'FontSize', fontsize);
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV');
%  title(sprintf('when coupling %.2f', coupling3))
% set(gca, 'linewidth',2, 'FontSize', fontsize);
%
% %%
%
% auc(1,:) = [auc_cg0  auc_np0  auc_sample0];
% auc(2,:) = [auc_cg  auc_np  auc_sample];
% auc(3,:) = [auc_cg2  auc_np2 auc_sample2];
% auc(4,:) = [auc_cg3  auc_np3 auc_sample3];
%
% % linespec = 'brk';
% %
% figure; bar(auc,'grouped');
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV');
% set(gca,'Ytick',[0:0.5:1]); set(gca,'YTickLabel',[0:0.5:1],'FontSize',fontsize);
% set(gca,'Xtick',[1:1:4]); set(gca,'XTickLabel',[0.07 0.1 0.2 0.3],'FontSize',fontsize);
%
% set(gca, 'linewidth',2, 'FontSize', fontsize);
% xlabel('\epsilon values','FontSize', fontsize)
% ylabel('Area under ROC','FontSize', fontsize)
%
% axis([0.5 4.5 0.5 1.01])
%
% figure;
% legend('Parametric PPLV','Non-parametric PPLV','Sample PLV');
% % text(1-.1,auc_cg+.03,sprintf('%.2f',auc_cg))
% % text(2-.1,auc_np+.03,sprintf('%.2f',auc_np))
% % text(3-.1,auc_sample+.03,sprintf('%.2f',auc_sample))













