%% Input

% fragMedian = [0.136 0.218 0.397];
% fragStd = [0.342 0.342 0.342];
% startDLRdata = [[1 20 70]/100, 0.1 0.1 0.1]';

fragMedian = [0.0576 0.1642 0.2301 0.3404];
fragStd = [0.3049 0.1924 0.2018 0.2041];

IMstripes = [0.052541 0.069785 0.096096 0.11577 0.1644 0.25819 0.29656 0.34869 0.40452 0.50184 0.58843];

load('LossGivenIM.mat')
totReconstructionCost = 1458000;
LRgivenIM = LossGivenIM / totReconstructionCost;

NsamplesLRgivenIM = 10000;
boundsDLR(:,1) = [0.05 0.25 0.6 0.8];
boundsDLR(:,2) = [0.15 0.4 0.73 1];
boundsCoV(:,1) = [0.05 0.05 0.05 0.05];
boundsCoV(:,2) = [1 1 1 1];

%% Run

calibrator = calibrateDLRs(fragMedian, fragStd, IMstripes, LRgivenIM, ...
    NsamplesLRgivenIM, boundsDLR, boundsCoV);

%calibrator = calibrator.nonLinearOptimisation;
calibrator = calibrator.particleSwarm;

%calibrator = calibrator.mapObjective(calibrator.DLRdataPSO);
%calibrator.plotObjective

calibrator.plotVulnerability
%calibrator.plotLRgivenIMdistributions('PSO')
