%% Input

fragMedian = [0.136 0.218 0.397];
fragStd = [0.342 0.342 0.342];

startDLRdata = [[10 20 30]/100, 0.3 0.3 0.3]';

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

startObj = objectiveFunction(startDLRdata, ...
                        IMstripes, LRgivenIM, empiricalMoments, ...
                        fragMedian, fragStd, ...
                        Nsamples, weightMoments);

%% Optimisation

% A*X <= B
A1 = -eye(numel(fragMedian)) + diag(ones(1, numel(fragMedian)-1), -1);
A2 = zeros(numel(fragMedian));
A3 = -eye(numel(fragMedian));

A = [A1, A2; A2, A3];
B = zeros(numel(startDLRdata),1);

% Aeq*X = Beq
Aeq = [];
Beq = [];

% bounds
lowerBound = zeros(numel(startDLRdata),1);
upperBound = [[0.2 0.7 1]'; 1*ones(numel(startDLRdata)/2,1)];

options = optimset('Display','iter');

objFun = @(DLRdata) objectiveFunction(DLRdata, ...
                        IMstripes, LRgivenIM, empiricalMoments, ...
                        fragMedian, fragStd, ...
                        Nsamples, weightMoments);

[finalDLRdata, finalObj, constraintsFlag] = fmincon(...
                        objFun, startDLRdata, ...
                        A, B, Aeq, Beq, lowerBound, upperBound, ...
                        @positiveAlphas, options);
        
finalDLRs = finalDLRdata(1:numel(startDLRdata)/2);
finalCoVdlrs = finalDLRdata(numel(startDLRdata)/2+1:end);

[~, CDFloss] = objectiveFunction([finalDLRs; finalCoVdlrs], ...
                        IMstripes, LRgivenIM, empiricalMoments, ...
                        fragMedian, fragStd, ...
                        Nsamples, weightMoments);

% Plot

cols = hsv(numel(IMstripes));

fragilities(:,1) = linspace(0, 2.5, 1000);
for ds = numel(fragMedian) : -1 : 1
    fragilities(:,ds+1) = logncdf(...
        fragilities(:,1), log(fragMedian(1,ds)), fragStd(1,ds));
end

finalVulnerability = VULNERABILITYbuilding(...
    fragilities, [0; finalDLRs], 'plot', 'S_{a}(T_{1}) [g]');

%figure; hold on
for im = 1 : numel(IMstripes)
    scatter(IMstripes(im)*ones(Nsamples,1), LRgivenIM(:,im), 5, 'k')
end
xlabel('IM [g]')
ylabel('Loss Ratio [-]')
set(gca, 'FontSize' ,18)


% figure; hold on
% for im = 1 : numel(IMstripes)
%     %histogram(LRgivenIM(:,im), 'normalization', 'pdf', 'DisplayStyle', 'stairs')
%     [empCDF(:,2), empCDF(:,1)] = ecdf(LRgivenIM(:,im));
%     plot(empCDF(:,1), empCDF(:,2), '--', 'Color', cols(im,:))
%     plot(CDFloss.LOSSdef, CDFloss.CDFlossIM(im,:), 'Color', cols(im,:))
%     
%     clear empCDF
% end
