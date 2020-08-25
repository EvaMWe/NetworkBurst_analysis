% this function i to select the short bursts that are not included in the
% Clusters or Agglomerates


function [shortStart, shortStop, nbInsiders] = selectShinLong(startShort, stopShort, startLong, stopLong)

Bstart = repmat(startShort',length(startLong),1); 
Log1 = Bstart >= startLong;
Log2 = Bstart <= stopLong;
Log = Log1.*Log2;

nbInsiders = sum(Log,2);

Log(sum(Log,2) == 0,:) = [];
Log2 = sum(Log,1);

shortStart = startShort(Log2~=0);
shortStop = stopShort(Log2~=0);
end
% get associated stop