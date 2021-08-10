function [c,ceq] = positiveAlphas(DLRdata,~,~,~,~,~,~,~)

Nds = numel(DLRdata)/2;
for ds = Nds : -1 : 1
    c(ds) = DLRdata(ds) + DLRdata(ds+Nds)^2 * DLRdata(ds) - 1;
end

ceq = [];

end