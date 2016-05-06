
clc; clear all; close all;

addpath(genpath(['C:\Users\Sergul\Documents\My Dropbox\multivariatePhase\matlabCodes']))

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes']))


MC = 1000;
couplingVec = [0.05 0.1 0.15 0.2 0.25];
sigmaVec =[0.5 1 1.5 2  2.5];
nVec =[1000 3000 5000 7000 9000];

% for idxCoupling = 1:length(couplingVec);
%     idxCoupling
    coupling = couplingVec(3);
%     for idxSigma = 1:length(sigmaVec);
%         idxSigma
        sigma = sigmaVec(3);
%         for idxN = 1:length(nVec);
%             idxN
            n = nVec(3);
            
            [FPcg TPcg FPsample TPsample auc_cg auc_sample] =...
                f_singleROCallTechnicalNote2(n,coupling, sigma, MC);

            
%         end
%     end
% end

save('singleROCtechnicalNote.mat','FPcg','TPcg','FPsample','TPsample','auc_cg','auc_sample','coupling','sigma','n');
