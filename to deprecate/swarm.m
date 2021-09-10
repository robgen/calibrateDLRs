Nvars = 2*numel(fragMedian);

objFun = @(DLRdata) objectiveFunction(DLRdata, ...
                        IMstripes, LRgivenIM, empiricalMoments, ...
                        fragMedian, fragStd, ...
                        Nsamples, weightMoments);

lowerBound = 0.001 * ones(numel(startDLRdata),1);
upperBound = [[0.2 0.7 0.99]'; 1*ones(numel(startDLRdata)/2,1)];

finalDLRdataPSO = particleswarm(objFun, Nvars, lowerBound, upperBound);
