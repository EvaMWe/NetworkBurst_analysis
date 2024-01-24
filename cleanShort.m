% this function i to select the short bursts that are not included in the
% Clusters or Agglomerates


function [shortStart, shortStop] = cleanShort(Startcluster, Stopcluster, burstStart, burstStop)

if length(Startcluster) <= 1 || length(Stopcluster) <= 1
    shortStart = burstStart;
    shortStop = burstStop;
    return
end

if length(Startcluster) ~= length(Stopcluster)
    [Startcluster, Stopcluster] = reSeq(Startcluster,Stopcluster);
end

if length(burstStart) ~= length(burstStop)
    [burstStart, burstStop] = reSeq(burstStart,burstStop);
end

if size(burstStart,2) > 1
    burstStart = burstStart';
end

if size(burstStop,2) > 1
    burstStop = burstStop';
end
Bstart = repmat(burstStart',length(Startcluster),1); 
Log1 = Bstart < Startcluster;
Log2 = Bstart > Stopcluster;
Log = Log1 + Log2;
sumLog = sum(Log,1);
shStartlog = sumLog == length(Startcluster);
shortStart = burstStart(shStartlog);
shortStop = burstStop(shStartlog);

if length(shortStart) ~= length(shortStop)
    [shortStart,shortStop] = reSeq(shortStart,shortStop);
end
end
% get associated stop