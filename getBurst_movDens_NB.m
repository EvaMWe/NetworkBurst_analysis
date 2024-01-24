% this is a new approach too improve the detection of fast bursts; longer
% periods without bursting should be ephasizied, it's more related to the
% real time line:
% use time bins and count the number of spikes within.

function [bStart_spkid, bStop_spkid,burstletStart, burstletStop] = getBurst_movDens_NB(spikeTime)

%
spS = 1;
stepWidth = 36;
sampleBin = 120;
spikeData = spikeTime*12000;
endSample = spikeData(end) -sampleBin;
firstSample = spikeData(1);
openWin = firstSample:stepWidth:endSample;
closeWin = openWin + sampleBin-1;
%WinMatrix = [openWin; closeWin];
%spikeTime_ms = spikeTime *1000;

%% count number of spikes per win
nbSpikes = arrayfun(@(a,b) countSpikes(a,b,spikeData), openWin,closeWin);


%% threshold crossing method: simply looking for the bins that cross a threshold of spikes as burst start and after returning below treshold dedicated as stops
logSpk = nbSpikes > spS;
logSpkDiff = diff(logSpk);
bStart = find(logSpkDiff == 1) +1;
bStop = find(logSpkDiff == -1);


bStart_t = openWin(bStart)./12000;
bStop_t = closeWin(bStop)./12000;

burstletStart = arrayfun(@(a) (find(spikeTime>a,1,'first')),bStart_t);
burstletStop= arrayfun(@(a) (find(spikeTime<a,1,'last')),bStop_t);

if burstletStop(1) < burstletStart(1)
    burstletStop = burstletStop(2:end);
end

if length(burstletStop) < length(burstletStart)
    burstletStart = burstletStart(1:length(burstletStart));
end

%% get clusters: this method goes along with the derivative; fast increasing sense (spike numbers in bins) are starters, 
fnc = (smooth(nbSpikes))';
SE = [1 1 1];
grad = [-2 -1 0 1 2];

fnc_grad1 = imfilter(fnc,grad);

%remove NaNs
if sum(isnan(fnc_grad1)) ~= 0
    idx = find(isnan(fnc_grad1));
    A = fnc_grad1(idx-1);
    B = fnc_grad1(idx +1);
    replacer = (B-A)/2+A;
    fnc_grad1(idx) = replacer;
end

dil = imdilate(fnc_grad1,SE);
startBurst = (dil - fnc_grad1) == 0;

h = mean(fnc_grad1) + std(fnc_grad1);
starters = fnc_grad1 >= 0.2*h & startBurst;

fnc_grad1inv = -fnc_grad1;
dil = imdilate(fnc_grad1inv,SE);
stopBurst = (dil - fnc_grad1inv) == 0;

h = mean(fnc_grad1inv) + std(fnc_grad1inv);
stops = fnc_grad1inv >= 0.2*h & stopBurst;

idxStart = find(starters == 1);
idxStop = find(stops == 1);
%% check for proper stops: condition: 15 following bin with a spk number < 15;

% check = arrayfun(@(follower) valCheck(follower,nbSpikes), idxStop);
% idxStop = idxStop(check);
% if idxStart(end) > idxStop(end)
%     idxStop = [idxStop length(stops)];
% end
[idxStart, idxStop] = reSeq2(idxStart,idxStop,nbSpikes);

%% recalculate to spike ID (from numbered time bin)

bStart_t = openWin(idxStart)./12000;
bStop_t = closeWin(idxStop)./12000;

bStart_spkid = arrayfun(@(a) (find(spikeTime>a,1,'first')),bStart_t);
bStop_spkid= arrayfun(@(a) (find(spikeTime<a,1,'last')),bStop_t);
end


function nbSpikes = countSpikes(a,b,spikeTime)
 nbSpikes = sum(spikeTime >= a & spikeTime < b);
end

function valid = valCheck(a, nbspikes)
follower = a:a+15;
if length(nbspikes) < follower(end)
    valid = false;
    return
end
valid = sum(nbspikes(follower) <= 1) == length(follower); %if smaller than length(a) at least one bin above threshold of 1 
end



% merge if needed
%S = Maximum_median(abs(filtered)(1:end-length(mask)+2),40,'Type','percent','Dimension',2);
%IntraBurst = burstStart(1,2:end)-burstStop(1,1:end-1);

%% visualize
%visualizeBurst(rateDiff, burstStart, burstStop, spikeTime);







