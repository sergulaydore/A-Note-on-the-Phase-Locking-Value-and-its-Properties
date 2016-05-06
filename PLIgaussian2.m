clc; clear all; close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes']))

bVec=[0.01:0.01:10];
lengthall = 100000;
for idxb = 1 %:length(bVec)

b = bVec(idxb);
b = 10
A = [1 b; b 1];

% Simulations

n = randn(2,500); % white Gaussian noise
n = n - repmat(mean(n,2),1,500);
x = A*n; % data generation
z = hilbert([x(1,:)' x(2,:)']);
A*A'*2
abs(cov([z(:,1) z(:,2)]))
Rzmat =inv(cov([z(:,1) z(:,2)]))
phase_mat_orig = angle(z);
diff_phase = phase_mat_orig(:,1)-phase_mat_orig(:,2);
PLI(idxb) = abs(mean(sign(diff_phase)));
PLV(idxb) = abs(mean(exp(j*diff_phase)));
k12(idxb) = -Rzmat(1,2)/(sqrt(Rzmat(1,1)*Rzmat(2,2)));
    end
% figure;
% plot(k12,PLV,'*')
% figure;
% plot(k12,PLI,'*')

zmat = [z(:,1) z(:,2)];
pCovSergul = abs((transpose(zmat)*zmat))/lengthall;
CovSergul = ((zmat)'*zmat)/lengthall
% figure;
% fs = 15;
% imagesc(k12Vec,mu12Vec,PLI'); caxis([0 1]); colorbar
% xlabel('Concentration Parameter','Fontsize',fs);
%  ylabel('Mean Offset','Fontsize',fs)
% 
%  figure;
% fs = 15;
% imagesc(k12Vec,mu12Vec,PLV'); caxis([0 1]); colorbar
% xlabel('Concentration Parameter','Fontsize',fs);
%  ylabel('Mean Offset','Fontsize',fs)

