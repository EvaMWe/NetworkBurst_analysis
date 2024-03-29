% this function performs accoding to moving Burst Frequency:
% Inter burst intervals are calculated and moving BUrst Frequency (rate) is
% calcualted) --- algorithm accoding to getBurst_movFR_NB

function [lBurstStart,lBurstStop,shStart,shStop] = getBlocks_slope(burstStart,burstStop,spikeVec, varargin)
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
IBI_s = (smooth(IBI_s(1:end)))';
mask = [1 1 1];
SE = [1 1 1];

if cluster == 1
    grad = [-1 -1 0 1 1];
else
    grad = [-2 -1 0 1 2];
end

%% find rough burst positions
filtered = conv(IBI_s,mask,'same');

rate = length(mask)./filtered;
rateDiff = imfilter(rate,grad);
%rateDiff = imfilter(rate,grad);


dil = imdilate(rateDiff,SE);
startBurst = (dil - rateDiff) == 0;
%startBurst = (dil -rateDiff) == 0;
%k�nnte noch durch einen noise filter ersetzt werden // statistik -->
%butterworth high pass, mean + 3*std = threshold;
filterSettings  = constructFilt;
[~,cleaner] =  filterSignal(rateDiff,filterSettings,'high' );

starters = rateDiff >= alpha_start*cleaner & startBurst;
%starters = rateDiff >= alpha_start*cleaner & startBurst;

rateDiffinv = -rateDiff;
dil = imdilate(rateDiffinv,SE);
stopBurst = (dil - rateDiffinv) == 0;
stops = rateDiffinv >= alpha_stop*cleaner & stopBurst;


%
idxStart_ = find(starters == 1);
idxStop_ = find(stops == 1);

%% refine starts // refines where developed with the rate curve...so with filtered it has to be the other way round

if isempty(idxStart_) || isempty(idxStop_)
    lBurstStart = [];
    lBurstStop = [];
    shStart = burstStart;
    shStop = burstStop;
    return
end

idxStart = refineStarts(idxStart_,rate);
idxStop = refineStops(idxStop_,rate);

if cluster == 0    
    [idxStart, idxStop]=reSeq2(idxStart',idxStop',rate);
else
    [idxStart, idxStop]=reSeq(idxStart',idxStop');
end

lBurstStart_ = burstStart(idxStart);
lBurstStop_ = burstStop(idxStop);

[lBurstStart, lBurstStop] = broadening_beta(lBurstStart_, lBurstStop_, spikeVec);

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
    Bstart = repmat(burstStart',length(lBurstStart),1);
    Log1 = Bstart < lBurstStart;
    Log2 = Bstart > lBurstStop;
    Log = Log1 + Log2;
    sumLog = sum(Log,1);
    shStartlog = sumLog == length(lBurstStart);
    shStart = burstStart(shStartlog);
    % get associated stop
    shStop = burstStop (shStartlog);
end





