% this function performs according to moving Burst Frequency:
% Inter burst intervals are calculated and moving BUrst Frequency (rate) is
% calcualted) --- algorithm accoding to getBurst_movFR_NB

function [lBurstStart,lBurstStop,shStart,shStop] = getBurstCluster_IBI(burstStart,burstStop,spikeVec, varargin)
%default
alpha_start = 1.5;
alpha_stop = 0.5;


if nargin > 4
    alpha_start = varargin{1};
    alpha_stop = varargin{2};
end



if length(burstStart) ~= length(burstStop)
    [burstStart,burstStop] = reSeq(burstStart,burstStop);
end

IBI_s = spikeVec(burstStart(2:end)) - spikeVec(burstStop(1:end-1));
%IBI_s = (smooth(IBI_s))';
mask = [1 1 1];

grad = [-2 0 2];


%% find rough burst positions
filtered = conv(IBI_s,mask,'same');
nBin = floor(length(burstStart)/5);
[nIBI, edges] = histcounts(filtered,nBin);
%[nIBI, edgesISI] = histcounts(IBI_s,350);
%figure,histogram('BinEdges',edges,'BinCounts',nrate)
idxCMA = calcCMA(nIBI);
% idxCMA_IBI = calcCMA(nIBI);
% h_IBI = edgesISI(6*(idxCMA_IBI+1));
h_min = edges(2*idxCMA+1);
h_max = edges(2*idxCMA+3);
filtered_smth = smooth(filtered,3)';
filteredDiff = imfilter(filtered_smth ,grad);

%% get burst starts = reflected as turning points from high values to low
% invert °1 derivative
filteredDiff_inv = -filteredDiff;
SE = [1 1 1];
filteredDiff_inv(isnan(filteredDiff_inv) | isinf(filteredDiff_inv)) = 0;
dil = imdilate(filteredDiff_inv,SE);

startBurst = (dil -filteredDiff_inv) == 0;
cleaner = std(filteredDiff_inv);

starters = filteredDiff_inv >= alpha_start*cleaner & startBurst;

%% get burst stops = reflected as turning points from low values to high
% °1 derivative
dil = imdilate(filteredDiff,SE);
stopBurst = (dil - filteredDiff) == 0;
stops = filteredDiff >= alpha_stop*cleaner & stopBurst;
%
idxStart_ = find(starters == 1)+1;
idxStop_ = find(stops == 1)-1;

%% refine starts
if isempty(idxStart_) || isempty(idxStop_)
    lBurstStart = [];
    lBurstStop = [];
    shStart = burstStart;
    shStop = burstStop;
    return
end

idxStart = unique(refineRight(idxStart_,filtered, h_max, h_min));
idxStop = unique(refineLeft(idxStop_,filtered, h_max, h_min));
if idxStop(end) > length(filtered)
    idxStop(end) = length(filtered);
end

if idxStart(end) > length(filtered)
    idxStart(end) = length(filtered);
end

if idxStop(1) == 0
    idxStop(1) = 1;
end

if idxStart(1) == 0
    idxStart(1) = 1;
end

%% check vality
[idxStart, idxStop]=reSeq(idxStart',idxStop');

logi = arrayfun(@(a) checkingStart(a, filtered, h_max), idxStart);
idxStart = idxStart(logi);

logi = arrayfun(@(a) checkingStop(a, filtered, h_max), idxStop);
idxStop = idxStop(logi);

[idxStart, idxStop]=reSeq(idxStart',idxStop');



logi = arrayfun(@(a,b) checkingBurst(a,b, filtered, h_max), idxStart, idxStop);
logi = logi == 0;
idxStop = idxStop(logi);
idxStart = idxStart(logi);

%     idxStart = arrayfun(@(a) shiftStart(IBI_s,a,h_IBI), idxStart);
%     idxStop = arrayfun(@(a) shiftStop(rate,a,h_max), idxStop);







%% refine starts according to IBI


lBurstStart = burstStart(idxStart)';
lBurstStop = burstStop(idxStop)';

%[lBurstStart, lBurstStop] = broadening_beta(lBurstStart_, lBurstStop_, spikeVec);

if length(lBurstStart) < 2 || length(lBurstStop) < 2
    shStart = burstStart;
    shStop = burstStop;
    return
end

if length(lBurstStart) ~= length(lBurstStop) 
    [lBurstStart,lBurstStop] = reSeq(lBurstStart,lBurstStop);
end

    
    %% select short bursts
    % starts

    Bstart = repmat(burstStart,length(lBurstStart),1);
    Log1 = Bstart < lBurstStart;
    Log2 = Bstart > lBurstStop;
    Log = Log1 + Log2;
    sumLog = sum(Log,1);
    shStartlog = sumLog == length(lBurstStart);
    shStart = burstStart(shStartlog);
    % get associated stop
    shStop = burstStop (shStartlog);
end

function valStart = checkingStart(idx,rate, maxi)
valStart = rate(idx) <= maxi;
end

function valStop = checkingStop(idx,rate, maxi)
valStop = rate(idx) <= maxi;
end

function exclude = checkingBurst(start,stop, rate, maxi)
exclude = rate(1,start:stop) >= maxi;
exclude = sum(exclude)/(stop-start);
exclude = exclude > 0.2;
end




