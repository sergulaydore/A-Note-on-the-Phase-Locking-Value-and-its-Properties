
% Here, we compare the estimated phase locking values with the analytical
% true value in Gaussian simulations

clc;
clear all;
close all;

addpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes\GaussianSimulations'])

N = 1000; % Monte Carlo

% Frequency range
fstart = 8;
delta_f = 2;

order = 50; % order of the filter
Fs = 100; % Sampling frequency
lengthall = 1000; % length of the signals
counter = 0; % counting signals which results NaN

length2vec =[5:5:100];
length(length2vec)
bVec = [0.128 0.34 0.65 0.899]; %0.193; %0.7325 %

% PLVsample = ones(length(bVec),length(length2vec),N);
% PLVsampleUb = ones(length(bVec),length(length2vec),N);
% PLVvonmises = ones(length(bVec),length(length2vec),N);
% PLVcircgauss = ones(length(bVec),length(length2vec),N);
% PLVanalCircgauss= ones(length(bVec),length(length2vec),N);
% plv = ones(length(bVec),length(length2vec),N);
% 
% fend = fstart + delta_f;
% freqrange = [fstart fend];
% h = filter_design('bandpass',freqrange,order,Fs,0);
% 
% % A = [1 0.8185; 0.8185 1];
% 
% for idxB = 1:length(bVec)
%     idxB
%     b = bVec(idxB);
%     A = [1 b; b 1];
%     
%     % Analytical Correlation results of X
%     Rmat_anal = A*A';
%     Rx1x2anal = Rmat_anal(1,2);
%     Rx1x1anal = Rmat_anal(1,1);
%     Rx2x2anal = Rmat_anal(2,2);
%     
%     % Analytical Correlation results of Y
%     hh = conv(h,fliplr(h));
%     hhh = conv(hh, fliplr(hh)); % due to filtfilt function
%     convresult = conv(Rx1x2anal,hhh);
%     Ry1y2anal = convresult(101);
%     convresult = conv(Rx1x1anal,hhh);
%     Ry1y1anal = convresult(101);
%     convresult = conv(Rx2x2anal,hhh);
%     Ry2y2anal = convresult(101) ;
%     
%     % Analytical Correlation results of Z (complex vector)
%     Rz1z1anal = 2*Ry1y1anal;
%     Rz2z2anal = 2*Ry2y2anal;
%     Rz1z2anal = 2*Ry1y2anal;
%     
%     Kmat_anal = [Rz1z1anal Rz1z2anal; Rz1z2anal Rz2z2anal];
%     invKmat_anal = inv(Kmat_anal);
%     
%     k12 = abs(invKmat_anal(1,2));
%     k11 = abs(invKmat_anal(1,1));
%     k22 = abs(invKmat_anal(2,2));
%     w = k12^2/(2*k11*k22-k12^2);
%     PLVcg = pi/sqrt(2)/k12^2*abs(det(invKmat_anal))*(w^(3/2)*hypergeom([3/4,5/4],1,w^2)+ 3/4*w^(5/2)*hypergeom([5/4,7/4],2,w^2))
%     
%     for mcindex = 1:N
%         mcindex;
%         %      A = [1 0.7*randn;0.7*randn 1]; % mixing matrix
%         
%         z = randn(2,lengthall); % white Gaussian noise
%         z = z - repmat(mean(z,2),1,lengthall);
%         x = A*z; % data generation
%         
%         % Empirical Correlation results of X
%         Rmat_emp = cov(x(1,:)',x(2,:)');
%         Rx1x2emp = Rmat_emp(1,2);
%         Rx1x1emp = Rmat_emp(1,1);
%         Rx2x2emp = Rmat_emp(2,2);
%         
%         % Filtered signal Y
%         
%         y(1,:) = filter_apply(x(1,:),h); %conv(xorig,h); %
%         y(2,:) = filter_apply(x(2,:),h); %conv(yorig,h); %
%         
%         
%         % Empirical Correlation results of Z (complex vector)
%         [invKmat_emp, r1_orig, r2_orig, rvec_orig, phase_mat_orig, z_orig] = f_amp_phase(y(1,:),y(2,:));
%         
%         while (isnan(det(invKmat_emp))),
%             
%             z = randn(2,lengthall); % white Gaussian noise
%             z = z - repmat(mean(z,2),1,lengthall);
%             x = A*z; % data generation
%             
%             % Filtered signal Y
%             fend = fstart + delta_f;
%             freqrange = [fstart fend];
%             h = filter_design('bandpass',freqrange,order,Fs,0);
%             y(1,:) = filter_apply(x(1,:),h); %conv(xorig,h); %
%             y(2,:) = filter_apply(x(2,:),h); %conv(yorig,h); %
%             [invKmat_emp, r1_orig, r2_orig, rvec_orig, phase_mat_orig, z_orig] = f_amp_phase(y(1,:),y(2,:));
%         end
%         
%         
%         % Empirical Correlation results of Y
%         covmat2 = cov(y(1,:),y(2,:));
%         Ry1y2emp = covmat2(1,2);
%         Ry1y1emp = covmat2(1,1);
%         Ry2y2emp = covmat2(2,2);
%         
%         
%         for lengthindex = 1:length(length2vec)
%             
%             length2 = length2vec(lengthindex);
%             y_part = zeros(2,length2);
%             PLVanalCircgauss(idxB,lengthindex,mcindex) = abs(PLVcg);
%             
%             % Estimate PLV over some samples
%             y_part(1,:) = y(1,lengthall/10:lengthall/10+length2-1); %+ 0.02*rand(1,length2+1);
%             y_part(2,:) = y(2,lengthall/10:lengthall/10+length2-1); %+ 0.02*rand(1,length2+1);
%             
%             % Circular Gaussian
%             [invKmat_par, r1, r2, rvecorg, phase_mat, z] = f_amp_phase(y_part(1,:),y_part(2,:)); % parameters of the distribution in polar coordinates
%             
%             k12 = abs(invKmat_par(1,2));
%             k11 = abs(invKmat_par(1,1));
%             k22 = abs(invKmat_par(2,2));
%             D = k12^2/(k11*k22);
%             w = k12^2/(2*k11*k22-k12^2);
%             PLVcg_par =  pi/(sqrt(2))*(1-1/D)* ...
%                 [ w^(3/2)*hypergeom([3/4 5/4],1, w^2) + ...
%                 3/4* w^(5/2)* hypergeom([5/4 7/4],2, w^2) ];
%             PLVcircgauss(idxB,lengthindex,mcindex) = abs(PLVcg_par);
%             
%             % Estimate PLV over all samples
%             [diff_phase] = f_diff_phase(phase_mat);
%             PLVsample(idxB,lengthindex,mcindex) = abs(mean(exp(sqrt(-1)*diff_phase)));
%             PLVsampleUb(idxB,lengthindex,mcindex) =...
%                 sqrt((PLVsample(idxB,lengthindex,mcindex)^2*length2-1)/(length2-1));
%             
%             %         % Von Mises
%             %         Kvm = f_fitmodel_cadieu(phase_mat); % Apply Cadieu's method to estimate coupling
%             %         coupling_sm(mcindex,lengthindex) = abs(Kvm(1,2));
%             coupling_ml = circ_kappa(phase_mat(:,1)-phase_mat(:,2)); % abs(Kvm(1,2)); %
%             PLVvonmises(idxB,lengthindex,mcindex) = abs(besseli(1,coupling_ml)/besseli(0,coupling_ml));
%             %         plv_par_vmsm(mcindex,lengthindex) = abs(besseli(1,coupling_sm(mcindex,lengthindex))/besseli(0,coupling_sm(mcindex,lengthindex)));
%             
%         end
%         
%     end
% end
% 
% save PLVsample PLVsample
% save PLVsampleUb PLVsampleUb
% save PLVvonmises PLVvonmises
% save PLVcircgauss PLVcircgauss
% save PLVanalCircgauss PLVanalCircgauss 

load PLVsample
load PLVsampleUb
load PLVvonmises
load PLVcircgauss
load PLVanalCircgauss

meanPLVsample = mean(PLVsample,3);
meanPLVsampleUb = mean(PLVsampleUb,3);
meanPLVvonmises = mean(PLVvonmises,3);
meanPLVcircgauss = mean(PLVcircgauss,3);
meanPLVanalCircgauss = mean(PLVanalCircgauss,3);

for idxB = 1:length(bVec);
    for idxLength = 1:length(length2vec)
        varPLVsample(idxB, idxLength) =...
            var(squeeze(PLVsample(idxB, idxLength,:)));
        
        varPLVsampleUb(idxB, idxLength) =...
            var(squeeze(PLVsampleUb(idxB, idxLength,:)));
        
        varPLVvonmises(idxB, idxLength) =...
            var(squeeze(PLVvonmises(idxB, idxLength,:)));
        
        varPLVcircgauss(idxB, idxLength) =...
            var(squeeze(PLVcircgauss(idxB, idxLength,:)));
    end
end

figure;
lwVec = [1.5 2.5 3 4.5]; fs = 15; ms = 8;
for idxB = [1 3]
    plot(length2vec,meanPLVsampleUb(idxB,:),'b-d','LineWidth',lwVec(idxB),'MarkerSize',ms);
hold on;
 plot(length2vec,meanPLVsample(idxB,:),'r-s','LineWidth',lwVec(idxB));
%plot(length2vec,meanPLVvonmises(idxB,:),'r-*','LineWidth',lw,'MarkerSize',ms);
plot(length2vec,meanPLVcircgauss(idxB,:),'k-o','LineWidth',lwVec(idxB),'MarkerSize',ms);
% plot(length2vec,meanPLVanalCircgauss(idxB,:),'m--','LineWidth',lw)
xlabel('Number of Observations','Fontsize',fs)
ylabel('Mean','Fontsize',fs)
% legend('Sample','Sample UB','Von Mises','Circ Gauss')

 h = legend('$\mbox{E} \left[ \sqrt{ \hat{\mbox{{PLV}}}^2_{\mbox{ub\_sample}} } \right]$',...
      '$\mbox{E}\left[\hat{\mbox{PLV}}_{\mbox{sample}}\right]$',...
     '$\mbox{E}\left[\hat{\mbox{PLV}}_{\mbox{circgauss}}\right]$');
 set(h,'Interpreter','latex','fontsize',13)
set(h,'box','off')
axis([0 100 0 1.6])
end

figure;
lw = 2; fs = 15;  ms = 8;
for idxB = [ 1 3]
    plot(length2vec,log10(varPLVsampleUb(idxB,:)),'b-d','LineWidth',lwVec(idxB));
hold on;
 plot(length2vec,log10(varPLVsample(idxB,:)),'r-s','LineWidth',lwVec(idxB));
%plot(length2vec,log10(varPLVvonmises(idxB,:)),'r.-','LineWidth',lw);
plot(length2vec,log10(varPLVcircgauss(idxB,:)),'k-o','LineWidth',lwVec(idxB));
% plot(length2vec,meanPLVanalCircgauss(idxB,:),'m--','LineWidth',lw)
xlabel('Number of Observations','Fontsize',fs)
ylabel('Variance (log)','Fontsize',fs)
 h = legend('$\log_{10} \left(\mbox{Var} \left[ \sqrt{ \hat{\mbox{PLV}}^2_{\mbox{ub\_sample}} } \right] \right)$',...
     '$\log_{10} \left( \mbox{Var}\left[\hat{\mbox{PLV}}_{\mbox{sample}}\right] \right)$',...
     '$\log_{10} \left( \mbox{Var}\left[\hat{\mbox{PLV}}_{\mbox{circgauss}}\right] \right)$');
 set(h,'Interpreter','latex','fontsize',13)
 set(h,'box','off')
axis([0 100 -3.5 0])
end