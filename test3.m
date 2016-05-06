clc;
clear all;
close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes']))

nVec =[5:5:100];
bVec =  [0.128 0.34 0.65 0.899]; %0.193; %0.7325 %
MC = 10000;

for idxB = 1:length(bVec)
    b = bVec(idxB);
    A = [1 b; b 1];
    
    % Analytical Correlation results of X
    Rmat_anal = A*A';
    Ry1y2anal = Rmat_anal(1,2);
    Ry1y1anal = Rmat_anal(1,1);
    Ry2y2anal = Rmat_anal(2,2);
    
    % Analytical Correlation results of Y
    %     hh = conv(h,fliplr(h));
    %     hhh = conv(hh, fliplr(hh)); % due to filtfilt function
    %     convresult = conv(Rx1x2anal,hhh);
    %     Ry1y2anal =convresult(101);
    %     convresult = %conv(Rx1x1anal,hhh);
    %     Ry1y1anal =convresult(101);
    %     convresult = %conv(Rx2x2anal,hhh);
    %     Ry2y2anal =convresult(101) ;
    
    % Analytical Correlation results of Z (complex vector)
    Rz1z1anal = 2*Ry1y1anal;
    Rz2z2anal = 2*Ry2y2anal;
    Rz1z2anal = 2*Ry1y2anal;
    
    Kmat_anal = [Rz1z1anal Rz1z2anal; Rz1z2anal Rz2z2anal];
    invKmat_anal = inv(Kmat_anal);
    
    k12 = abs(invKmat_anal(1,2));
    k11 = abs(invKmat_anal(1,1));
    k22 = abs(invKmat_anal(2,2));
    w = k12^2/(2*k11*k22-k12^2);
    PLVcg = pi/sqrt(2)/k12^2*abs(det(invKmat_anal))*(w^(3/2)*hypergeom([3/4,5/4],1,w^2)+ 3/4*w^(5/2)*hypergeom([5/4,7/4],2,w^2))
    
    for idxN = 1:length(nVec)
        idxN
        Length = nVec(idxN);
        
        for idxMC = 1:MC
            z = randn(2,Length); % white Gaussian noise
            z = z - repmat(mean(z,2),1,Length);
            x = A*z; % data generation
            
            % Empirical Correlation results of X
            Rmat_emp = cov(x(1,:)',x(2,:)');
            Rx1x2emp = Rmat_emp(1,2);
            Rx1x1emp = Rmat_emp(1,1);
            Rx2x2emp = Rmat_emp(2,2);
            
            [invKmat_emp, r1_orig, r2_orig, rvec_orig, phase_mat, z_orig] =...
                f_amp_phase(x(1,:),x(2,:));
            
            k12 = abs(invKmat_emp(1,2));
            k11 = abs(invKmat_emp(1,1));
            k22 = abs(invKmat_emp(2,2));
            D = k12^2/(k11*k22);
            w = k12^2/(2*k11*k22-k12^2);
            PLVcg_par =  pi/(sqrt(2))*(1-1/D)* ...
                [ w^(3/2)*hypergeom([3/4 5/4],1, w^2) + ...
                3/4* w^(5/2)* hypergeom([5/4 7/4],2, w^2) ];
            PLVcircgauss(idxB,idxN,idxMC) = abs(PLVcg_par);
            
             [diff_phase] = f_diff_phase(phase_mat);
            PLVsample(idxB,idxN,idxMC) = abs(mean(exp(sqrt(-1)*diff_phase)));
            PLVsampleUb(idxB,idxN,idxMC) =...
               sqrt((PLVsample(idxB,idxN,idxMC)^2*Length-1)/(Length-1));
           
             coupling_ml = circ_kappa(phase_mat(:,1)-phase_mat(:,2)); % abs(Kvm(1,2)); %
            PLVvonmises(idxB,idxN,idxMC) =...
                abs(besseli(1,coupling_ml)/besseli(0,coupling_ml));
        end
    end
end


meanPLVsample = mean(PLVsample,3);
meanPLVsampleUb = mean(PLVsampleUb,3);
meanPLVvonmises = mean(PLVvonmises,3);
meanPLVcircgauss = mean(PLVcircgauss,3);
% meanPLVanalCircgauss = mean(PLVanalCircgauss,3);

for idxB = 1:length(bVec);
    for idxN = 1:length(nVec)
        varPLVsample(idxB, idxN) =...
            var(squeeze(PLVsample(idxB, idxN,:)));
        
        varPLVsampleUb(idxB, idxN) =...
            var(squeeze(PLVsampleUb(idxB, idxN,:)));
        
        varPLVvonmises(idxB, idxN) =...
            var(squeeze(PLVvonmises(idxB, idxN,:)));
        
        varPLVcircgauss(idxB, idxN) =...
            var(squeeze(PLVcircgauss(idxB, idxN,:)));
    end
end

save PLVsample PLVsample
save PLVsampleUb PLVsampleUb
save PLVvonmises PLVvonmises
save PLVcircgauss PLVcircgauss
% save PLVanalCircgauss PLVanalCircgauss 

% load PLVsample
% load PLVsampleUb
% load PLVvonmises
% load PLVcircgauss

figure;
lw = 3; fs = 15;
for idxB = 1:length(bVec)
    plot(nVec,meanPLVsampleUb(idxB,:),'b--','LineWidth',lw);
hold on;
plot(nVec,meanPLVsample(idxB,:),'r--','LineWidth',lw);
plot(nVec,meanPLVvonmises(idxB,:),'r','LineWidth',lw);
plot(nVec,meanPLVcircgauss(idxB,:),'k--','LineWidth',lw);
% plot(length2vec,meanPLVanalCircgauss(idxB,:),'m--','LineWidth',lw)
xlabel('Number of Observations','Fontsize',fs)
legend('Sample','Sample UB','Von Mises','Circ Gauss')

h = legend('$\mbox{E} \left[ \sqrt{ {PLV}^2_{ub\_sample} } \right]$',...
    '$\mbox{E}\left[{PLV}_{sample}\right]$', '$\mbox{E}\left[{PLV}_{vonmises}\right]$',...
    '$\mbox{E}\left[{PLV}_{circgauss}\right]$');
set(h,'Interpreter','latex')

% axis([0 100 -.2 1.2])
end

figure;
lw = 3; fs = 15;
for idxB = 1:length(bVec)
    plot(nVec,log10(varPLVsampleUb(idxB,:)),'b','LineWidth',lw);
hold on;
plot(nVec,log10(varPLVsample(idxB,:)),'r--','LineWidth',lw);
plot(nVec,log10(varPLVvonmises(idxB,:)),'r','LineWidth',lw);
plot(nVec,log10(varPLVcircgauss(idxB,:)),'k--','LineWidth',lw);
% plot(length2vec,meanPLVanalCircgauss(idxB,:),'m--','LineWidth',lw)
xlabel('Number of Observations','Fontsize',fs)
h = legend('$\log_{10} \left(\mbox{Var} \left[ \sqrt{ {PLV}^2_{ub\_sample} } \right] \right)$',...
    '$\log_{10} \left(\mbox{Var}\left[{PLV}_{sample}\right] \right)$', '$\log_{10} \left( \mbox{Var}\left[{PLV}_{vonmises}\right] \right)$',...
    '$\log_{10} \left( \mbox{Var}\left[{PLV}_{circgauss}\right] \right)$');
set(h,'Interpreter','latex')
end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            