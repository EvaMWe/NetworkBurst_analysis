%This function is used to reanaylse the spike train regarding
%Burstagglomerates: "fast consecutive bursts including the very fast
%consecutive bursts designated as burstclsuter"

function [StartAggl, StopAggl, restStart, restStop] = findAgglomerats(shBurstStart,shBurstStop, Startcluster, Stopcluster, spikeVec)

allStarters = sort([shBurstStart; Startcluster]);
allStops = sort([shBurstStop; Stopcluster]);

%allStarters = LBstart;
%allStops = LBstop;

[StartAggl,StopAggl] = getBlocks_slope(allStarters, allStops,spikeVec,0.3, 0.3);

if length(StartAggl) ~= length(StopAggl)
    [StartAggl,StopAggl] = reSeq(StartAggl,StopAggl);
end

[StartAggl, StopAggl] = mergeBursts(StartAggl, StopAggl, spikeVec, 1, 10);
[restStart,restStop] = cleanShort(StartAggl, StopAggl, allStarters, allStops);

if length(restStart) ~= length(restStop)
    [restStart,restStop] = reSeq(restStart,restStop);
end



end