clc; clear all; close all;

addpath(genpath('F:\29.02.2012backup\Research Spring 2012\February\ISBI_2012\matlab codes\Sergul'));
addpath(genpath('F:\29.02.2012backup\Research Spring 2012\april\paper_draft\matlab_codes'))
load('F:\29.02.2012backup\Research Spring 2012\useful_old_stuff\dropbox\group_AR\MonkeyData\data\DataAllGo'); % 1 x 18 cell

group1 = [1 2 3 4 10 13 14 17];
group2 = [5 6 7 8 9 11 12 15 16 18];
idxTime = 42;
channelsInteresting = [3 5];

freqLow = 13;
freqHigh = 22;

idxTimeVec = [idxTime-5:1:idxTime+5]; % 120 msec
dataAll = [];
for k = 1:18 %length(group1)
    idxExp = k
    loadTimeFreqDataCompare
    
    nFreq = size(timeFreqGo,2);
    nTime= size(timeFreqGo,3);
    
    nTrialGo(k) = size(timeFreqGo,4);
    nTrialNoGo(k) = size(timeFreqNoGo,4);
    
%     for m = 1:nTime
%         for idxFreq = 1:nFreq
%             data = squeeze(timeFreqGo(:,idxFreq,m,:));
%             nTrial = size(data,2);
%             meanData = mean(data,2);
%             stdData = std(data')';
%             dataZeroMean = (data-repmat(meanData,1,nTrial));
%             dataProc = dataZeroMean./repmat(stdData,1,nTrial);
%             dataAll = [dataAll dataProc];
%         end
%     end
    
end

% phase_mat = angle(dataAll);
% diff_phase_mat12 = squeeze(phase_mat(1,:)-phase_mat(2,:));
% diff_phase_mat12(diff_phase_mat12<-pi) = diff_phase_mat12(diff_phase_mat12<-pi) + 2*pi;
% diff_phase_mat12(diff_phase_mat12>pi) = diff_phase_mat12(diff_phase_mat12>pi) - 2*pi;
% 
% x = [-pi:0.01:pi];
% invKmat12 = inv(cov([(dataAll(1,:))' (dataAll(2,:))']));
% y = f_pdf_marginal(invKmat12,x);
% 
% coupling_ml = circ_kappa(diff_phase_mat12); % abs(Kvm(1,2)); %
% mean_ml = circ_mean(diff_phase_mat12);
% y_vm = 1/(2*pi*besseli(0,coupling_ml))*exp(coupling_ml * cos(x - mean_ml));
% 
% 
% N = size(phase_mat,2);
% for index = 1:N
%     index
%     plv_cg(index) = randarb(x,y);
%     plv_vm(index) = randarb(x,y_vm);
%     
% end
% figure
% qqplot_sergul(diff_phase_mat12,plv_cg,[0 0 1]); % title(sprintf('%d and %d',cx,cy)); %xlabel('x bok'); ylabel('y bok')
% hold on;
% plot([-3:0.01:3],[-3:0.01:3],'--r','LineWidth', 2.0)