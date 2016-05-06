
% clc; clear all; close all;

dataDir =  ['F:\29.02.2012backup\Research Spring 2012\april\'...
        'paper_draft\matlab_codes\ResultsSignificance\lowBeta\biPLV120\']; %\PPLVschelter
    
    % parameters
alpha = 0.1;
idx = 1;
nChannels = 15;
thr =3;
disp_x = 7;
disp_y =15;

PPLVgoCount = zeros(nChannels,nChannels);
PPLVnogoCount = zeros(nChannels,nChannels);

% expVec =  [1 2 3 4 10 13 14 17]; %
%  expVec =  [5 6 7 8 9 11 12 15 16 18];
expVec  =[1:1:18];
% expVec = [2 5 6 8 13 14 18]
valuesAll = [];
for idxExp =expVec
     nstr = num2str(idxExp);
    
    if idxExp<10
        nstr = ['0' nstr];
    end
    
    filename = [dataDir 'sigTest' nstr '.mat'];
    load(filename);

    [IndX IndY] = find(ones(15,15)-(1-triu(ones(15,15)))-eye(15)==1);
    for pIdx = 1:length(IndX)
        valuesGo(pIdx) = PPLVgoActual(IndX(pIdx),IndY(pIdx));
        valuesNoGo(pIdx) = PPLVnogoActual(IndX(pIdx),IndY(pIdx));
    end
    
    valuesAll = [valuesAll valuesGo valuesNoGo];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataDir =  ['F:\29.02.2012backup\Research Spring 2012\april\'...
        'paper_draft\matlab_codes\ResultsSignificance\lowBeta\biPLV260\']; %\PPLVschelter
    
    % parameters
alpha = 0.1;
idx = 1;
nChannels = 15;
thr =3;
disp_x = 7;
disp_y =15;

PPLVgoCount = zeros(nChannels,nChannels);
PPLVnogoCount = zeros(nChannels,nChannels);

% expVec =  [1 2 3 4 10 13 14 17]; %
%  expVec =  [5 6 7 8 9 11 12 15 16 18];
expVec  =[1:1:18];
% expVec = [2 5 6 8 13 14 18]
% valuesAll = [];
for idxExp =expVec
     nstr = num2str(idxExp);
    
    if idxExp<10
        nstr = ['0' nstr];
    end
    
    filename = [dataDir 'sigTest' nstr '.mat'];
    load(filename);

    [IndX IndY] = find(ones(15,15)-(1-triu(ones(15,15)))-eye(15)==1);
    for pIdx = 1:length(IndX)
        valuesGo(pIdx) = PPLVgoActual(IndX(pIdx),IndY(pIdx));
        valuesNoGo(pIdx) = PPLVnogoActual(IndX(pIdx),IndY(pIdx));
    end
    
    valuesAll = [valuesAll valuesGo valuesNoGo];
    
end