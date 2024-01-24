% this function performs accoding to moving Burst Frequency:
% Inter burst intervals are calculated and moving BUrst Frequency (rate) is
% calcualted) --- algorithm accoding to getBurst_movFR_NB

function [lBurstStart,lBurstStop,shStart,shStop] = getBurstCluster_slope(burstStart,burstStop,spikeVec, varargin)
%default
alpha_start = 1.5;
alpha_stop = 0.5;
cluster = 0; %special adjustments for clustesr detection



if nargin > 4
    alpha_start = varargin{1};
    alpha_stop = varargin{2};
end

if nargin == 6
    cluster = varargin{3};
end

if length(burstStart) ~= length(burstStop)
    [burstStart,burstStop] = reSeq(burstStart,burstStop);
end

IBI_s = spikeVec(burstStart(2:end)) - spikeVec(burstStop(1:end-1));
%IBI_s = (smooth(IBI_s))';
mask = [1 1 1 1 1];

grad = [-2 0 2];


%% find rough burst positions
filtered = conv(IBI_s,mask,'same');

rate_ = length(mask)./filtered;
f = length(burstStart)/200;
if f > 1
     %% refine baseline
    [~,td_baseline_] = tempBase(rate_, f, 0 );
    td_baseline = baseRef(td_baseline_, [1 1 1 1], [-2 0 2]);
    rate_0 = (rate_ - td_baseline)./td_baseline;
    rate = rate_0;
    rate(rate < 0) = 0;
    [nrate, edges] = histcounts(rate,200);
    %[nIBI, edgesISI] = histcounts(IBI_s,200);
    %figure,histogram('BinEdges',edges,'BinCounts',nrate)
    idxCMA = calcCMA(nrate(2:end));
   % idxCMA_IBI = calcCMA(nIBI);
   % h_IBI = edgesISI(6*(idxCMA_IBI+1));
   if 6*(idxCMA+1) > 0.5*length(edges)
       h_max = edges(round(0.5*length(edges)));
   else
       h_max = edges(6*(idxCMA+1));
   end
   
    factor = td_baseline./Minimum_median(td_baseline, 50, 'Type', 'percent', 'Dimension', 2);
    hDyn = factor.*h_max;
    %% calculate a dynamic threshold
    
    rateSm = smooth(rate)';
    rateDiff = imfilter(rateSm,grad);

    SE = [1 1 1 1];
else
    rate = rate_;    
    hDyn = repmat(Maximum_median(rate, 60, 'Type', 'percent', 'Dimension', 2),1,length(rate));
    %h_IBI = Minimum_median(IBI_s, 40, 'Type', 'percent', 'Dimension', 2);
    grad = [-2 0 2];
    rateDiff = imfilter(rate,grad);
    rateDiff(1) = mean(rateDiff(2:3));
    rateDiff(end) = mean(rateDiff(end-3:end-2));
    SE = [1 1 1];
    cluster = 1;
end

rateDiff(isnan(rateDiff) | isinf(rateDiff)) = 0;
dil = imdilate(rateDiff,SE);

startBurst = (dil -rateDiff) == 0;
cleaner = std(rateDiff);

starters = rateDiff >= alpha_start*cleaner & startBurst;

rateDiffinv = -rateDiff;
dil = imdilate(rateDiffinv,SE);
stopBurst = (dil - rateDiffinv) == 0;
stops = rateDiffinv >= alpha_stop*cleaner & stopBurst;
%
idxStart_ = find(starters == 1);
idxStop_ = find(stops == 1);

%% refine starts
if isempty(idxStart_) || isempty(idxStop_)
    lBurstStart = [];
    lBurstStop = [];
    shStart = burstStart;
    shStop = burstStop;
    return
end

if cluster == 0
    idxStart = unique(refineStarts(idxStart_,rate, hDyn,1));
    idxStop = unique(refineStops(idxStop_,rate,hDyn,1));
    if idxStop(end) > length(rate)
        idxStop(end) = length(rate);
    end
    
    %% check vality
    if idxStart(1) == 0
        idxStart(1) =1;
    end
    logi = arrayfun(@(a) checkingStart(a, rate, hDyn), idxStart);
    idxStart = idxStart(logi);
    
    logi = arrayfun(@(a) checkingStop(a, rate, hDyn), idxStop);
    idxStop = idxStop(logi);
    
%     idxStart = arrayfun(@(a) shiftStart(IBI_s,a,h_IBI), idxStart);
%     idxStop = arrayfun(@(a) shiftStop(rate,a,h_max), idxStop);
    [idxStart, idxStop]=reSeq(idxStart',idxStop');
    
else
    idxStart = idxStart_;
    idxStop = idxStop_;
    idxStop = arrayfun(@(a) shiftStop(rate,a, hDyn), idxStop);
    [idxStart, idxStop]=reSeq(idxStart',idxStop');
    
end



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
valStart = rate(idx) <= maxi(idx);
end

function valStop = checkingStop(idx,rate, maxi)
valStop = rate(idx) <= maxi(idx);
end

function starts = shiftStart(inter,xStart,h)
x = xStart;
while inter(x) > h(xStart)
    x = x+1;
end
starts = x;
end

function stop = shiftStop(inter,xStop,h)
x = xStop;
while inter(x) > h(xStop)
    x = x+1;
end
stop = x;
end




