%% Input

fragMedian = [0.136 0.218 0.397];
fragStd = [0.342 0.342 0.342];
startDLRdata = [[1 20 70]/100, 0.1 0.1 0.1]';

IMstripes = [0.034424 0.070101 0.10882 0.14801 0.20314 0.24688 0.28296 0.33767 0.43135 0.55804];

load('LossGivenIM.mat')
totReconstructionCost = 1458000;
LRgivenIM = LossGivenIM / totReconstructionCost;

%% Run

calibrator = calibrateDLRs(fragMedian, fragStd, IMstripes, LRgivenIM);

calibrator = calibrator.nonLinearOptimisation;
calibrator = calibrator.particleSwarm;

calibrator = calibrator.mapObjective(calibrator.DLRdataPSO);
calibrator.plotObjective

calibrator.plotVulnerability
calibrator.plotLRgivenIMdistributions('PSO')
