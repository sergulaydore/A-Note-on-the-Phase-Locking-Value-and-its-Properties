clc; clear all; close all;

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes']))

k12Vec = [0.01:0.1:0.9];
mu12Vec = [0:0.001:pi];
lengthall = 500;
for idxk12 = 1:length(k12Vec)
    for idxmu12 = 1:length(mu12Vec)

k12 = k12Vec(idxk12);
mu12 = mu12Vec(idxmu12);
Korig = [1 -k12*exp(j*mu12) ; -k12*exp(-j*mu12) 1]/(1-k12^2);
Aorig = chol((Korig));

% Simulations

n = randn(2,lengthall); % white Gaussian noise
n = n - repmat(mean(n,2),1,lengthall);
x = Aorig*n; % data generation
z = hilbert([x(1,:)' x(2,:)']);
Rzmat =inv(cov([z(:,1) z(:,2)]));
phase_mat_orig = angle(z);
diff_phase = phase_mat_orig(:,1)-phase_mat_orig(:,2);
PLI(idxk12,idxmu12) = abs(mean(sign(diff_phase)));
PLV(idxk12,idxmu12) = abs(mean(exp(j*diff_phase)));
    end
end


figure;
fs = 15;
imagesc(k12Vec,mu12Vec,PLI'); caxis([0 1]); colorbar
xlabel('Concentration Parameter','Fontsize',fs);
 ylabel('Mean Offset','Fontsize',fs)

 figure;
fs = 15;
imagesc(k12Vec,mu12Vec,PLV'); caxis([0 1]); colorbar
xlabel('Concentration Parameter','Fontsize',fs);
 ylabel('Mean Offset','Fontsize',fs)

