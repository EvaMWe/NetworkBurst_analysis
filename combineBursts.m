% select bursts from Var(-> getBurstCluster_slope) that are not inside
% bursts from Fixed (-> getBurstCluster_IBI)
%
% in contrast to selectShinLong: all bursts coming from
% getBurstCLuster_slope that starts or ends within bursts coming from
% getBurstCLuster_IBI are discarded

function [burstStarts, burstStops] = combineBursts(startFixed, stopFixed, startVar, stopVar) 
if size(startVar,2) == 1
    startVar = startVar';
    stopVar =stopVar';
end

if size(startFixed,2) ~= 1
    startFixed = startFixed';
    stopFixed =stopFixed';
end

size_2=max(stopVar(end), stopFixed(end));
LogiMat = zeros(length(startVar), size_2);
%Create burst matrix; one burst per row
for line = 1:length(startVar)
    LogiMat(line,startVar(line):stopVar(line)) = 1;
end

LogoIdx = arrayfun(@(startFixed, stopFixed) getLog(startFixed, stopFixed,LogiMat), startFixed, stopFixed);


%%put together
LogoId = LogoIdx;
LogoId(LogoId == 0) = [];
Template = zeros(1,length(startVar));
Template(LogoId) = 1;
burstIdx = find(Template == 0);
burstStarts = sort([startFixed; startVar(burstIdx)']);
burstStops = sort([stopFixed; stopVar(burstIdx)']);
end


function idx = getLog(start, stop, M)
logVec = zeros(1,size(M,2));
logVec(1, start:stop) = 1;
logMat = repmat(logVec, size(M,1), 1);
idxLog = logMat == 1 & M == 1;
idx = sum(idxLog,2);
idx = find(idx~=0);
if isempty(idx)
    idx=0;
elseif size(idx,1) ~= 1
    idx = idx(1);
end
end