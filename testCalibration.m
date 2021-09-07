%% Generate dummy LRgivenIM data based on a set of DLRs

fragMedian = [0.1 0.2 0.4];
fragStd = [0.3 0.3 0.3];
targetDLRdata = [[3 20 80]/100, 0.1 0.5 0.3]';

IMstripes = linspace(0.05,0.6,5);

Nsamples = 1000;

% irrelevant - only needed to run the objective function once
LRgivenIM = rand(Nsamples,numel(IMstripes));
Nsamples = size(LRgivenIM, 1);
empiricalMoments = rand(4, numel(IMstripes));
% irrelevant - only needed to run the objective function once

weightMoments = ones(size(empiricalMoments,1),1);

% derive LRgivenIM distribution
[~, inputCDFloss] = objectiveFunction(targetDLRdata, ...
                        IMstripes, LRgivenIM, empiricalMoments, ...
                        fragMedian, fragStd, ...
                        Nsamples, weightMoments);

% simulate LR consistent with the input LRgivenIM distribution
rng(1)
unifRand = rand(Nsamples,numel(IMstripes));
for im = 1 : numel(IMstripes)
    LRgivenIM(:,im) = interp1([0 inputCDFloss.CDFlossIM(im,:)], ...
        [0; eps; inputCDFloss.LOSSdef(2:end)], unifRand(:,im));
end
   
for im = numel(IMstripes) : -1 : 1
    empiricalMoments(1,im) = mean(LRgivenIM(:,im));
    empiricalMoments(2,im) = var(LRgivenIM(:,im));
    empiricalMoments(3,im) = skewness(LRgivenIM(:,im));
    empiricalMoments(4,im) = kurtosis(LRgivenIM(:,im));
end

%% Run optimisation based on the dummy data

close all
startDLRdata = rand(6,1);
optimiseDLRs
plotDLRs

%% Verify 