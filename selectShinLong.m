% this function i to select the short bursts that are not included in the
% Clusters or Agglomerates


function [shortStart, shortStop, nbInsiders] = selectShinLong(startShort, stopShort, startLong, stopLong)
%%%CHECK FOR DIMENSIONS!!!
if size(startShort,2) ==1
    startShort = startShort';
end
if size(stopShort,2) ==1
    stopShort = stopShort';
end

Bstart = repmat(startShort,length(startLong),1); 
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