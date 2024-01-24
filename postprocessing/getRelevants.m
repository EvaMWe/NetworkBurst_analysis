function relF = getRelevants(coeff,eigVal,varNames)
%% make scee calculation to get number of principle components
EV = cumsum(eigVal);
eigValrel = EV./EV(end);
idx = find(eigValrel>0.8);
nbPC = idx(1);
varStorage = cell(nbPC*length(varNames),1);

count = 1;
for c = 1:nbPC
    %get contribution of individual features to PC
    coeff_ = coeff(:,c);
    Csq = coeff_.^2;
    port = sort(Csq, 'descend');
    
    
    %gest most relevant features
    portCum = cumsum(port);
    idx = find(portCum>0.8);
    nbPort = idx(1);
    
    % isolate parameters
    names = varNames(1:nbPort);
    varStorage(count:count+length(names)-1) = names;
    count = count + length(names);
end

idxCell = ~(cellfun(@isempty, varStorage));
relF = varStorage(idxCell);
end
