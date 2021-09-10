% final vulnearbility
fragilities(:,1) = linspace(0, 2.5, 1000);
for ds = numel(fragMedian) : -1 : 1
    fragilities(:,ds+1) = logncdf(...
        fragilities(:,1), log(fragMedian(1,ds)), fragStd(1,ds));
end
finalVulnerability = VULNERABILITYbuilding(...
    fragilities, [0; finalDLRs], 'NOplot', 'S_{a}(T_{1}) [g]');


figure('Position', [280   27   560   420]); hold on
plot(finalVulnerability(:,1), finalVulnerability(:,2), 'LineWidth', 2)
scatter(fragMedian, ...
    interp1(finalVulnerability(:,1), finalVulnerability(:,2), fragMedian), ...
    100, 'filled')

for im = 1 : numel(IMstripes)
    scat = scatter(IMstripes(im)*ones(Nsamples,1), LRgivenIM(:,im), 5, 'k');
end
set(gca, 'XLim', [0 1.5*max(IMstripes)])
xlabel('IM [g]')
ylabel('Loss Ratio [-]')
set(gca, 'FontSize' ,18)

for ds = numel(fragMedian) : -1 : 1
    strDLR{ds} = sprintf(...
        'DLR_{%d}=%1.2f%%(%1.2f)', ds, finalDLRs(ds), finalCoVdlrs(ds));
end
annotation(gcf,'textbox',...
    [0.54742857142857 0.140476190476191 0.313285714285716 0.288095238095239],...
    'String',strDLR,'LineStyle','none',...
    'FontSize',18,'FontName','Helvetica Neue','FitBoxToText','off');


im = randi(numel(IMstripes));
clear empCDF
figure('Position', [841   27   560   420]); hold on
%histogram(LRgivenIM(:,im), 'normalization', 'pdf', 'DisplayStyle', 'stairs')
[empCDF(:,2), empCDF(:,1)] = ecdf(LRgivenIM(:,im));
plot(empCDF(:,1), empCDF(:,2), '--', 'Color', 'k', 'LineWidth', 2)
plot(CDFloss.LOSSdef, CDFloss.CDFlossIM(im,:), 'Color', 'k', 'LineWidth', 2)
xlabel('Loss Ratio, lr [-]')
ylabel('P(LR\leqlr) [-]')
set(gca, 'FontSize' ,18)
title(sprintf('IM=%2.4f', IMstripes(im)))