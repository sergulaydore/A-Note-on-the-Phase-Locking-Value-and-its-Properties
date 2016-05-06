
% clc;
% clear all;
% close all;
% n = 10000;
% MC = 5000;

function [TPcg12 FPcg TPnp12 FPnp TPsample12 FPsample auc_cg auc_np...
    auc_sample] = f_ROCtechnicalNote(coupling, MC, n,sigma)

M = 2;
a = .15; b = .2; c=10;
w = [1.03 1.01];
  h = 0.02;
u0 = ones(3*M,1)/3*M;

e = zeros(M,M);

e(1,2) = coupling; e(2,1) = e(1,2);

delta_thresh = 0.01;
min_thresh = 0.0001; max_thresh = 1-0.0001; 
thresh_vec = [min_thresh:delta_thresh:max_thresh];
  

for mcindex = 1:MC
    mcindex
    noise = randn(ceil(n*4/3),3*M); noise(:,2:3:end) = 0; noise(:,3:3:end) = 0;
    x = bst_roessler(w, e, n, 1/h, a, b, c, sigma*noise(:,1:3:end), u0)';
    
    [invKmat_par, r1, r2, rvecorg, phase_mat, z] = f_amp_phase(x(1,:),x(2,:)); % parameters of the distribution in polar coordinates
            
            k12 = abs(invKmat_par(1,2));
            k11 = abs(invKmat_par(1,1));
            k22 = abs(invKmat_par(2,2));
            D = k12^2/(k11*k22);
            w = k12^2/(2*k11*k22-k12^2);
            PLVcg_par =  pi/(sqrt(2))*(1-1/D)* ...
                [ w^(3/2)*hypergeom([3/4 5/4],1, w^2) + ...
                3/4* w^(5/2)* hypergeom([5/4 7/4],2, w^2) ];
            PLVcircgauss(mcindex) = abs(PLVcg_par);
            
            % Estimate PLV over all samples
            [diff_phase] = f_diff_phase(phase_mat);
            PLVsample(mcindex) = abs(mean(exp(sqrt(-1)*diff_phase)));
            PLVsampleUb(mcindex) =...
                sqrt((PLVsample(mcindex)^2*n-1)/(n-1));
 
end

for idx_thresh = 1:length(thresh_vec)
    idx_thresh
    threshold = thresh_vec(idx_thresh);
 FPcg(idx_thresh) = length(find(PLVcg_23>threshold))/MC;
TPcg12(idx_thresh) = length(find(PLVcg_12>threshold))/MC;
% TPcg13(idx_thresh) = length(find(PLVcg_13>threshold))/MC;

FPnp(idx_thresh) = length(find(R23_1>threshold))/MC;
TPnp12(idx_thresh) = length(find(R12_3>threshold))/MC;
% TPnp13(idx_thresh) = length(find(R13_2>threshold))/MC;

FPsample(idx_thresh) = length(find(R23>threshold))/MC;
TPsample12(idx_thresh) = length(find(R12>threshold))/MC;
% TPsample13(idx_thresh) = length(find(R13>threshold))/MC;

end

auc_cg = trapz(sort(FPcg),sort(TPcg12));
auc_np = trapz(sort(FPnp),sort(TPnp12));
auc_sample = trapz(sort(FPsample),sort(TPsample12));

% save PLVcg_12 PLVcg_12;
% save PLVcg_23 PLVcg_23;
% 
% save R12_3 R12_3;
% save R23_1 R23_1;
% 
% save R12 R12;
% save R23 R23;

% fid = fopen('cg.txt', 'w');
% fprintf(fid,'Circular Gaussian \r\n Large \r\n');
% fprintf(fid,'%s\r\n ',PLVcg_23);
% fprintf(fid,'*');
% fprintf(fid,'%s\r\n ',PLVcg_12);
% fprintf(fid,'*');
% fclose(fid);



% figure; plot(log10(FPcg),(TPcg12),'*'); hold on; plot(log10(FPnp),(TPnp12),'*r');



% X = [ R23' PLVcg_23' R23_1'];
% p = [25 50 75];
% percentiles = prctile(X,p)'
% 
% % G ={{'Sample PLV'},{'Parametric PPLV'},{'Non-parametric PPLV'}}
% figure; linewidth =2.5; fontsize=20; % boxplot(X,G,'notch', 'on'); 
% boxplot(X,'notch', 'on'); 
%  a = findall(gca,'type','line'); % get all lines in the box plot
% for j=1:length(a);
%     set(a(j), 'linewidth',linewidth)
% end
% %  set(gca,'Xtick',[1:1:3]); set(gca,'XTickLabel',G);
% %  set(gca,'fontsize',fontsize);
% 
%  axis([0.5 3.5 0 0.7])
%  set(gca,'Ytick',[0:0.1:0.7]); set(gca,'YTickLabel',[0:0.1:0.7]);
%  set(gca,'fontsize',fontsize);
%  
%  set(findobj(gca,'Type','text'),'fontunits','points')
%  set(findobj(gca,'Type','text'),'fontsize',fontsize-7)
%  
% [intercept slope error]=f_test_mvn(x); % axis([0 40 0 40])
%  set(gca,'fontsize',fontsize);