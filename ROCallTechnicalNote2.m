clc; clear all; close all;

addpath(genpath(['C:\Users\Sergul\Documents\My Dropbox\multivariatePhase\matlabCodes']))

addpath(genpath(['F:\29.02.2012backup\Research Spring 2012\useful_old_stuff'...
    '\dropbox\Manuscript-2\matlab_codes']))


MC = 1000;

for idxMC = 1:MC
    idxMC
    M = 2;
    a = .15; b = .2; c=10;
    wRoessler = [1.03 1.01];
    sigma = 1.5;  h = 0.02;
    u0 = ones(3*M,1)/3*M;
    n = 10000;
    e = zeros(M,M);
    coupling = 0.15;
    
    KmatCoupled = NaN;
    KmatNull = NaN;
    
    while (isnan(det(KmatCoupled))) || (isnan(det(KmatNull)))
    e(1,2) = coupling; e(2,1) = e(1,2);
    noise = randn(ceil(n*4/3),3*M); noise(:,2:3:end) = 0; noise(:,3:3:end) = 0;
    xCoupled = bst_roessler(wRoessler, e, n, 1/h, a, b, c, sigma*noise(:,1:3:end), u0)';
    zCoupled = xCoupled - repmat(mean(xCoupled,2),1,n);
    [KmatCoupled, r1org, r2org, rvecorg, phaseMatCoupled, z] = ...
        f_amp_phase(zCoupled(1,:),zCoupled(2,:));
    
    k12 = abs(KmatCoupled(1,2));
    k11 = abs(KmatCoupled(1,1));
    k22 = abs(KmatCoupled(2,2));
    D = k12^2/(k11*k22);
    w = k12^2/(2*k11*k22-k12^2);
    PLVcg_par =  pi/(sqrt(2))*(1-1/D)* ...
        [ w^(3/2)*hypergeom([3/4 5/4],1, w^2) + ...
        3/4* w^(5/2)* hypergeom([5/4 7/4],2, w^2) ];
    PLVcgCoupled(idxMC) = abs(PLVcg_par);
    % Estimate PLV over all samples
    [diff_phase] = f_diff_phase(phaseMatCoupled);
    PLVsampleCoupled(idxMC) = abs(mean(exp(sqrt(-1)*diff_phase)));
    
    
    eNull = zeros(M,M);
    xNull = bst_roessler(wRoessler, eNull, n, 1/h, a, b, c, sigma*noise(:,1:3:end), u0)';
    zNull = xNull - repmat(mean(xNull,2),1,n);
    [KmatNull, r1org, r2org, rvecorg, phaseMatNull, z] = ...
        f_amp_phase(zNull(1,:),zNull(2,:));
    
    k12 = abs(KmatNull(1,2));
    k11 = abs(KmatNull(1,1));
    k22 = abs(KmatNull(2,2));
    D = k12^2/(k11*k22);
    w = k12^2/(2*k11*k22-k12^2);
    PLVcg_par =  pi/(sqrt(2))*(1-1/D)* ...
        [ w^(3/2)*hypergeom([3/4 5/4],1, w^2) + ...
        3/4* w^(5/2)* hypergeom([5/4 7/4],2, w^2) ];
    PLVcgNull(idxMC) = abs(PLVcg_par);
    % Estimate PLV over all samples
    [diff_phase] = f_diff_phase(phaseMatNull);
    PLVsampleNull(idxMC) = abs(mean(exp(sqrt(-1)*diff_phase)));
     end
end

delta_thresh = 0.01;
min_thresh = 0.0001; max_thresh = 1-0.0001;
thresh_vec = [min_thresh:delta_thresh:max_thresh];

for idx_thresh = 1:length(thresh_vec)
    idx_thresh
    threshold = thresh_vec(idx_thresh);
    FPcg(idx_thresh) = length(find(PLVcgNull>threshold))/MC;
    TPcg(idx_thresh) = length(find(PLVcgCoupled>threshold))/MC;
    FPsample(idx_thresh) = length(find(PLVsampleNull>threshold))/MC;
    TPsample(idx_thresh) = length(find(PLVsampleCoupled>threshold))/MC;
end
    
figure;
plot(FPsample, TPsample)
hold on;
plot(FPcg,TPcg,'r')
    
