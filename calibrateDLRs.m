classdef calibrateDLRs
        
    properties
        fragMedian
        fragStd
        IMstripes
        LRgivenIM
        
        objFun
        
        mappedObj
        mappedMLR
        mappedCOV
        
        DLRdataNLO
        constraintsFlagNLO
        CDFsLossNLO
        
        DLRdataPSO
        CDFsLossPSO
    end
    
    properties(Access = 'private')
        A
        B
        Aeq
        Beq
        lowerBound
        upperBound
    end
    
    methods
        
        function self = calibrateDLRs(fragMedian, fragStd, ...
                IMstripes, LRgivenIM, NsamplesLRgivenIM)
            
            if nargin < 5; NsamplesLRgivenIM = 10000; end
            
            self.fragMedian = fragMedian;
            self.fragStd = fragStd;
            self.IMstripes = IMstripes;
            self.LRgivenIM = LRgivenIM;
            
            self = self.setConstraints;
            
            % irrelevant - to be removed after changing objFx
            empiricalMoments = rand(4, numel(self.IMstripes));
            weightMoments = ones(size(empiricalMoments,1),1);
            % irrelevant - to be removed after changing objFx
            
            self.objFun = @(DLRdata) objectiveFunction(DLRdata, ...
                self.IMstripes, self.LRgivenIM, empiricalMoments, ...
                self.fragMedian, self.fragStd, ...
                NsamplesLRgivenIM, weightMoments);
        end
        
        %% Macro methods
        
        function self = mapObjective(self, baselineDLRdata, Nseed)
            
            if nargin < 3; Nseed = 100; end
            
            for DS = numel(self.fragMedian) : -1 : 1
                
                [self.mappedMLR(:,:,DS), self.mappedCOV(:,:,DS)] = meshgrid(...
                    linspace(self.lowerBound(DS), self.upperBound(DS), Nseed), ...
                    linspace(self.lowerBound(numel(self.fragMedian)+DS), ...
                    self.upperBound(numel(self.fragMedian)+DS), Nseed));
                
                seedDLRdata = baselineDLRdata;
                for row = Nseed : -1 : 1
                    for col = Nseed : -1 : 1
                        seedDLRdata([DS,numel(self.fragMedian)+DS]) = ...
                            [self.mappedMLR(row,col,DS) ...
                             self.mappedCOV(row,col,DS)];
                        
                        self.mappedObj(row, col, DS) = ...
                            self.objFun(seedDLRdata);
                    end
                end
            end
        end
                
        function self = nonLinearOptimisation(self, options, startDLRdata)
            
            if nargin < 3; startDLRdata = rand(2*numel(self.fragMedian),1); end
            if nargin < 2
                options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
                    'OptimalityTolerance', 1e-12, ...
                    'Algorithm','interior-point', 'PlotFcn',{@optimplotx,...
                    @optimplotfval,@optimplotfirstorderopt});
            end
            
            [self.DLRdataNLO, ~, self.constraintsFlagNLO] = fmincon(...
                self.objFun, startDLRdata, self.A, self.B, ...
                self.Aeq, self.Beq, self.lowerBound, self.upperBound, ...
                @positiveAlphas, options);
            
            [~, self.CDFsLossNLO] = self.objFun(self.DLRdataNLO);
        end
        
        function self = particleSwarm(self, options)
            
            if nargin < 2
                options = optimoptions('particleswarm', ...
                    'Display', 'iter', 'FunctionTolerance', 10e-10);
            end
            
            Nvars = 2*numel(self.fragMedian);
            
            DLRdata = particleswarm(self.objFun, Nvars, ...
                self.lowerBound, self.upperBound, options);
            
            self.DLRdataPSO = DLRdata(:);
            
            [~, self.CDFsLossPSO] = self.objFun(DLRdata(:));
        end
        
        function plotObjective(self)
            
            if isempty(self.mappedObj)
                disp('Run mapObjective first')
                return
            end
            
            for DS = 1 : numel(self.fragMedian)
                MLR = self.mappedMLR(:,:,DS);
                COV = self.mappedCOV(:,:,DS);
                OBJ = self.mappedObj(:,:,DS);
                [~, indMin] = min(OBJ(:));
                                
                figure
                surf(MLR, COV, OBJ); hold on
                scatter3(MLR(indMin), COV(indMin), OBJ(indMin), 150, 'k', 'filled')
                xlabel(sprintf('MLR(DS_%d)', DS))
                ylabel(sprintf('CoV(DS_%d)', DS))
                zlabel('Objective function')
                set(gca, 'FontSize', 18)
            end
            
        end
        
        function plotVulnerability(self)
            
            if ~isempty(self.DLRdataNLO)
                finalVulnerability(:,1) = self.IMstripes;
                finalVulnerability(:,2) = mean(1-self.CDFsLossNLO.CDFlossIM, 2);
                
                figure; hold on
                for im = 1 : numel(self.IMstripes)
                    s = scatter(self.IMstripes(im)*ones(size(self.LRgivenIM,1),1), ...
                        self.LRgivenIM(:,im), 5, 'k', 'filled');
                end
                
                p = plot(finalVulnerability(:,1), finalVulnerability(:,2), ...
                    'r', 'LineWidth', 2);
                
                f = scatter(self.fragMedian, ...
                    interp1(finalVulnerability(:,1), finalVulnerability(:,2), ...
                    self.fragMedian), 100, 'filled');
                
                legend([s,p,f], {'Empirical', 'Vulnerability', 'Frag. Medians'})
                
                for ds = numel(self.fragMedian) : -1 : 1
                    strDLR{ds} = sprintf(...
                        'DLR_{%d}=%1.2f%%(%1.2f)', ds, self.DLRdataNLO(ds), ...
                        self.DLRdataNLO(numel(self.fragMedian)+ds));
                end
                annotation(gcf,'textbox',...
                    [0.54742857142857 0.140476190476191 0.313285714285716 0.288095238095239],...
                    'String',strDLR,'LineStyle','none',...
                    'FontSize',18,'FontName','Helvetica Neue','FitBoxToText','off');
                
                xlabel('IM [g]')
                ylabel('Loss Ratio [-]')
                set(gca, 'FontSize' ,18)
                
                title('Nonlinear Optimisation')
            else
                warning('Nonlinear optimisation was not run. Can''t plot it')
            end
            
            if ~isempty(self.DLRdataPSO)
                
                finalVulnerability(:,1) = self.IMstripes;
                finalVulnerability(:,2) = mean(1-self.CDFsLossPSO.CDFlossIM, 2);
                
                figure; hold on
                for im = 1 : numel(self.IMstripes)
                    s = scatter(self.IMstripes(im)*ones(size(self.LRgivenIM,1),1), ...
                        self.LRgivenIM(:,im), 5, 'k', 'filled');
                end
                
                p = plot(finalVulnerability(:,1), finalVulnerability(:,2), ...
                    'r', 'LineWidth', 2);
                
                f = scatter(self.fragMedian, ...
                    interp1(finalVulnerability(:,1), finalVulnerability(:,2), ...
                    self.fragMedian), 100, 'filled');
                
                legend([s,p,f], {'Empirical', 'Vulnerability', 'Frag. Medians'})
                
                for ds = numel(self.fragMedian) : -1 : 1
                    strDLR{ds} = sprintf(...
                        'DLR_{%d}=%1.2f%%(%1.2f)', ds, self.DLRdataPSO(ds), ...
                        self.DLRdataPSO(numel(self.fragMedian)+ds));
                end
                annotation(gcf,'textbox',...
                    [0.54742857142857 0.140476190476191 0.313285714285716 0.288095238095239],...
                    'String',strDLR,'LineStyle','none',...
                    'FontSize',18,'FontName','Helvetica Neue','FitBoxToText','off');
                
                xlabel('IM [g]')
                ylabel('Loss Ratio [-]')
                set(gca, 'FontSize' ,18)
                
                title('Particle Swarm Optimisation')
            else
                warning('Particle swarm optimisation was not run. Can''t plot it')
            end
            
        end
        
        function plotLRgivenIMdistributions(self, optMethod)
            
            for im = 1 : numel(self.IMstripes)
                
                clear empCDF
                
                figure('Position', [841   27   560   420]); hold on
                [empCDF(:,2), empCDF(:,1)] = ecdf(self.LRgivenIM(:,im));
                plot(empCDF(:,1), empCDF(:,2), '--', 'Color', 'k', 'LineWidth', 2)
                plot(self.(['CDFsLoss' optMethod]).LOSSdef, ...
                    self.(['CDFsLoss' optMethod]).CDFlossIM(im,:), ...
                    'Color', 'k', 'LineWidth', 2)
                xlabel('Loss Ratio, lr [-]')
                ylabel('P(LR\leqlr) [-]')
                set(gca, 'FontSize' ,18)
                legend('Empirical', 'Fitted', 'Location', 'SouthEast')
                title(sprintf('IM=%2.4f', self.IMstripes(im)))
            end
        end
        
        %% Micro methods
        
        function self = setConstraints(self)
            
            Nds = numel(self.fragMedian);
            
            % A*X <= B
            A1 = -eye(Nds) + diag(ones(1, Nds-1), -1);
            A2 = zeros(Nds);
            A3 = -eye(Nds);
            
            self.A = [A1, A2; A2, A3];
            self.B = zeros(2*Nds,1);
            
            % Aeq*X = Beq
            self.Aeq = [];
            self.Beq = [];
            
            % bounds
            self.lowerBound = 0.001 * ones(2*Nds,1);
            self.upperBound = [[0.2 0.7 0.99]'; 1*ones(Nds,1)]; warning('make me independent of the number of DSs')
            
        end
        
    end
end
