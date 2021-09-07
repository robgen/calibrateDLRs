function [obj, CDFsLoss] = objectiveFunction(DLRdata, IMstripes, ...
                                 empiricalLossGivenIMsamples, ...
                                 empiricalLossGivenIMmoments, ...
                                 fragMedian, fragStd, ...
                                 Nsamples, weightMoments)

meanDLRs = [0 DLRdata(1:numel(DLRdata)/2)'];
covDLRs = [0 DLRdata(numel(DLRdata)/2+1:end)'];

fragilities = getFragilities(IMstripes, fragMedian, fragStd);
CDFsLoss = getCDRlossGivenIM(fragilities, meanDLRs, covDLRs);
lossGivenIMsamples = getLRgivenIMsamples(CDFsLoss, Nsamples);

% % check with discrete KS test
% pValues = runKStests(lossGivenIMsamples, empiricalLossGivenIMsamples);
% obj = combinePvalues(pValues);

% check with central moments
maxOrderMoments = size(empiricalLossGivenIMmoments,1);
guessedMoments = getLRgivenIMmoments(lossGivenIMsamples, maxOrderMoments);
obj = combineMoments(guessedMoments, empiricalLossGivenIMmoments, weightMoments);

end


function fragilities = getFragilities(IMstripes, fragMedian, fragStd)

fragilities(:,1) = IMstripes;
for ds = numel(fragMedian) : -1 : 1
    fragilities(:,ds+1) = logncdf(...
        IMstripes, log(fragMedian(1,ds)), fragStd(1,ds));
end

end


function CDFsLoss = getCDRlossGivenIM(fragilities, meanDLRs, covDLRs)

% calculate probability of having different damage states, given IM
probDS = [ones(size(fragilities,1),1) fragilities(:,2:end)] - ...
    [fragilities(:,2:end) zeros(size(fragilities,1),1)];
% NOTE: the first column contains the probability of DS0 (no damage)

% P(L<l|DSk): beta  distributions
alfa = (1-meanDLRs)./covDLRs.^2 - meanDLRs;
if any(alfa < 0)
    alfa(alfa<0) = 0.01;
end

beta = alfa .* (1-meanDLRs) ./ meanDLRs;

LOSSdef = linspace(0,1,50+1)';

for ds = numel(meanDLRs) : -1 : 2
    CDFlossDSmatrix(1,:,ds) = betacdf(LOSSdef, alfa(ds), beta(ds));
end
CDFlossDSmatrix(1,:,1) = 1; % P(L<l|DS0) = 1 for any l


%%% write it avoiding the loop
for im = size(fragilities,1) : -1 : 1
    CDFlossDSmatrix(im,:,:) = CDFlossDSmatrix(1,:,:);
end

for l = numel(LOSSdef):-1:1; probDS3d(:,l,:) = probDS; end
%%% write it avoiding the loop

% calculate P(L<l|IM)
CDFsLoss.LOSSdef = LOSSdef;
CDFsLoss.CDFlossIM = sum( CDFlossDSmatrix .* probDS3d, 3 );

% add P(L<l|Dsk) to the structure
CDFsLoss.CDFlossDS = squeeze(CDFlossDSmatrix(1,:,:));

end


function lossGivenIMsamples = getLRgivenIMsamples(CDFsLoss, Nsamples)
    
for im = size(CDFsLoss.CDFlossIM,1) : -1 : 1
    rawECDF = [ CDFsLoss.LOSSdef, CDFsLoss.CDFlossIM(im,:)' ];
    [~,in] = unique(rawECDF(:,2));
    ECDF = rawECDF(in,:);
    
    lossGivenIMsamples(:,im) = ECDFsampler(ECDF, rand(Nsamples,1));
end

end


function ECDFsamples = ECDFsampler(ECDF, uniformRandomNumbers)

ECDF(1,1) = 10*eps;
ECDF = [0 0; ECDF];

onesInCDF = find(ECDF(:,2) > 1-eps);
ECDF(onesInCDF(2:end),:) = [];

ECDFsamples = interp1(ECDF(:,2), ECDF(:,1), uniformRandomNumbers, 'linear');

end


function lossGivenIMmoments = getLRgivenIMmoments(lossGivenIMsamples, maxOrder)
    
for im = size(lossGivenIMsamples,2) : -1 : 1
    lossGivenIMmoments(1,im) = mean(lossGivenIMsamples(:,im));
    lossGivenIMmoments(2,im) = var(lossGivenIMsamples(:,im));
    lossGivenIMmoments(3,im) = skewness(lossGivenIMsamples(:,im));
    lossGivenIMmoments(4,im) = kurtosis(lossGivenIMsamples(:,im));
end

end


function obj = combineMoments(guessedMoments, targetMoments, weightMoments)

targetVector = weightMoments .* (guessedMoments - targetMoments);
obj = sum(sum( targetVector.^2 ));

end


function pValues = runKStests(lossGivenIMsamples, empiricalLossGivenIMsamples)

for im = size(lossGivenIMsamples,2) : -1 : 1
    [~, pValues(im)] = kstest2(...
        lossGivenIMsamples(:,im), empiricalLossGivenIMsamples(:,im));
end

pValues = 10^100 * pValues;
end


function obj = combinePvalues(pValues)
    
    % minimise the negative = maximise
    obj = sum(-pValues);
end
