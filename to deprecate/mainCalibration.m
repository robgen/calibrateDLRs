%% Input

% fragMedian = [0.0493 0.1364 0.2177 0.3975];
% fragStd = 0.342 * [1 1 1 1];
% startDLRdata = [[1 2 3 4]/100, 0.1 0.1 0.1 0.1]';

fragMedian = [0.136 0.218 0.397];
fragStd = [0.342 0.342 0.342];
startDLRdata = [[1 20 70]/100, 0.1 0.1 0.1]';

IMstripes = [0.034424 0.070101 0.10882 0.14801 0.20314 0.24688 0.28296 0.33767 0.43135 0.55804];

load('LossGivenIM.mat')
totReconstructionCost = 1458000;
LRgivenIM = LossGivenIM / totReconstructionCost;

Nsamples = size(LRgivenIM, 1);
for im = numel(IMstripes) : -1 : 1
    empiricalMoments(1,im) = mean(LRgivenIM(:,im));
    empiricalMoments(2,im) = var(LRgivenIM(:,im));
    empiricalMoments(3,im) = skewness(LRgivenIM(:,im));
    empiricalMoments(4,im) = kurtosis(LRgivenIM(:,im));
end
weightMoments = ones(size(empiricalMoments,1),1);

%% Optimisation

close all
optimiseDLRs
plotDLRs
