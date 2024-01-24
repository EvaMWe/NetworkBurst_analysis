% this function i to select the short bursts that are not included in the
% Clusters or Agglomerates


function [outsideSpikes] = isolate_oSpikes(startBurst,stopBurst, spikeList)
if size(startBurst,2)>1
    startBurst = startBurst';
    stopBurst = stopBurst';
end
indexVec = (1:length(spikeList));
Bstart = repmat(indexVec ,length(startBurst),1); 
Log1 = Bstart < startBurst;
Log2 = Bstart > stopBurst;
Log = Log1 + Log2;
sumLog = sum(Log,1);
index = sumLog == length(stopBurst);
outsideSpikes = spikeList(index);

end
% get associated stop