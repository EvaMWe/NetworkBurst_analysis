function [burstStart, burstStop] = broadening_beta(start_temp, stop_temp, spikeTimes)

%delta = stop_temp-start_temp;
%f = @computeS;

%matrix(:,3) = arrayfun(@(a,b)max(diff(spikeTimes(a:b))),start_temp,stop_temp);

%%elongate



% enough starts and stops
if length(start_temp)<= 1 || length(stop_temp) <=1
    burstStart = start_temp;
    burstStop = stop_temp;
    return
%same length?
elseif length(start_temp) ~= length(stop_temp)
    [start_temp,stop_temp] = reSeq (start_temp,stop_temp);
    if length(start_temp) ~= length(stop_temp)
        burstStart = start_temp;
        burstStop = stop_temp;
        return
    end
end

stop_temp(stop_temp < start_temp(1)) = [];

if sum(stop_temp - start_temp < 0) ~= 0
    [start_temp,stop_temp] = reSeq (start_temp,stop_temp);
    if sum(stop_temp - start_temp < 0) ~= 0
        burstStart = start_temp;
        burstStop = stop_temp;
        return
    end
end
    
winLen = 100;
win = zeros(1,winLen);
hood = 5;
spikeTimesLong = [win spikeTimes];
startLong = start_temp + winLen;
stopLong = stop_temp + winLen;
burstStart = startDirection(startLong, stopLong, spikeTimesLong, winLen,hood);

spikeTimesLong = [spikeTimes win];
startLong = start_temp;
stopLong = stop_temp;
burstStop = stopDirection(stopLong, startLong, spikeTimesLong, winLen, hood);

end

%% START
function start = startDirection(start_temp, stop_temp, spikeTimes, cut, win)

if size(spikeTimes,2) > 1
    spikeTimes = spikeTimes';
end

if size(start_temp,2) > 1
    start_temp = start_temp';
end

if size(stop_temp,2) > 1
    stop_temp = stop_temp';
end

nbSp = zeros(length(start_temp),1);
summe = 1;
adding = ones(length(start_temp),1);

winStart = start_temp - win;
winStop = start_temp -1;
alphapre = arrayfun(@(a,b) computeS(a,b,spikeTimes), start_temp, stop_temp)...
    ./arrayfun(@(a,b) computeS(a,b,spikeTimes), winStart+1, winStop+1);
count = 1;

while summe ~= 0 
    
        
    IBIin = arrayfun(@(a,b) computeS(a,b,spikeTimes), start_temp-count+1, stop_temp); % clac MW of ISI within the bursts
        
    IBIpre = arrayfun(@(a,b) computeS(a,b,spikeTimes), winStart-count+1, winStop-count+1);
    
    alpha = IBIin./IBIpre;
    temp = alpha < alphapre;  
   
    
    %temp(start_temp(2:end)-count <= stop_temp(1:end-1)) = 0;
    logic = winStart(2:end)-count <= stop_temp(1:end-1);
    temp(logical([1; logic])) = 0;
    
    if start_temp(1)-count <= 0
        temp(1,1) = 0;
    end
    
    out = spikeTimes(start_temp(2:end)-count+1)<=spikeTimes(start_temp(1:end-1));
    temp(out) = 0;
    
    adding = temp & adding; % save the zero!
    summe = sum(adding);
    nbSp = nbSp+adding;
   
    alphapre = alpha;
    count = count+1;
    
  
end
start = start_temp - nbSp;
start = start - cut;
end

%% STOP
function stop = stopDirection(stop_temp,start_temp, spikeTimes, cut, win)


if size(spikeTimes,2) > 1
    spikeTimes = spikeTimes';
end

if size(start_temp,2) > 1
    start_temp = start_temp';
end

if size(stop_temp,2) > 1
    stop_temp = stop_temp';
end


nbSp = zeros(length(stop_temp),1);
summe = 1;
adding = ones(length(stop_temp),1);

winStart = stop_temp + 1;
winStop = stop_temp + win;
alphapre = 1000;
count = 1;

while summe ~= 0
    IBIin = arrayfun(@(a,b) computeS(a,b,spikeTimes), start_temp, stop_temp+count-1); % calc MW of ISI within the bursts
        
    IBIpost = arrayfun(@(a,b) computeS(a,b,spikeTimes), winStart+count-1, winStop+count-1);
    
    alpha = IBIin./IBIpost;
    temp = alpha < alphapre; 
    
    temp(winStop(1:end-1)+count >= start_temp(2:end,1)) = 0;
    
    if stop_temp(end,1) + count >= length(spikeTimes-cut)
        temp(1,1) = 0;
    end
    
    adding = temp & adding; % save the zero!
    summe = sum(adding);
    nbSp = nbSp+adding;
   
    alphapre = alpha;
    count = count+1;
    
end
stop = stop_temp + nbSp;
end

function thresh = computeS(a,b,spikeTimes)
thresh = max(diff(spikeTimes(a:b)));
end


