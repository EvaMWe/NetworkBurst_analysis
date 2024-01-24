function [burstInfo, burstID, start, stop] = burstExploration(start,stop,spikeVec,spikeTime,varargin)
%default
%theta = 10;
minContr = 5;
minSpk = 3;
nb_burst = length(start);
burstInfo = repmat(struct('burstNumber',1),nb_burst,1);
burstID = zeros(nb_burst,1);
count = 1;

if nargin >= 5
    minSpk = varargin{1};
end

if nargin == 6
    minContr = varargin{2};
end

if length(start) <= 1 || length(stop) <=1
    for burst = 1:nb_burst
        burstInfo(count).nbcontribChannels = 0;
        count = count +1;
    end
    
    if nb_burst == 0
        burstInfo(1).nbcontribChannels = 0;
    end
    return
end
    
if length(start) ~= length(stop)
    [start,stop] = reSeq(start,stop);
end

nb_burst = length(start);

if iscell(spikeTime)
    spikeTime = conversion(spikeTime,1);
end

%% for each

for burst = 1:nb_burst
    if start(burst) <= 0
        start(burst) = 1;
    end
    start_temp = spikeVec(start(burst));
    
    stop_temp = spikeVec(stop(burst));
    Alin = find(spikeTime >= start_temp & spikeTime <= stop_temp);
    [Arow, Acol] = ind2sub(size(spikeTime),Alin);
    contributors = unique(Acol);
    nbContr = length(contributors);
    nbSpikes = stop(burst)-start(burst) +1;
    burstInfo(count).burstNumber=burst;
    if nbContr >= minContr && nbSpikes >= nbContr*minSpk                
        burstInfo(count).burstSubscripts = [Arow,Acol];
        burstInfo(count).contribChannels = contributors;
        burstInfo(count).nbcontribChannels = nbContr;
        burstInfo(count).burstID = burst;
        burstID(count) = burst;
        burstInfo(count).nbSpikes = nbSpikes;
        count = count +1;
    else
        burstInfo(count).nbcontribChannels = 0;
        count = count +1;
    end
end

burstID(burstID == 0) = [];
start = start(burstID);
stop = stop(burstID);

if length(start) ~= length(stop)
    [start,stop] = reSeq(start,stop);
end

end