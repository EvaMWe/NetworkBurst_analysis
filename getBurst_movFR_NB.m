%core algorithm for network burst detection.
% input: spikeTime, spiketimes as a vector
% (1) find rough burst positions:
%     moving firing rate: convolution of spike Time with mask --> time
%     within a window defined by a fixed number of spikes (length(mask)
%     rate: number of spikes / time 
%     rateDiff: 1. derivative of rate --> strong increases indicates start
%     of burst, strong decreases indicates stop of a burst
%     filter rateDiff to get noise and threshold to clean data -->
%     idxStart,idxStop
%  
% (2) rearrange idxStart and idxStop

function [burstStart, burstStop] = getBurst_movFR_NB(spikeTime)


%defaults
mask = [1 0 0 0 0 -1];
SE = [1 1 1 1 1];
grad = [-2 -1 0 1 2];

%% find rough burst positions
filtered = conv(spikeTime,mask,'same');
if sum(filtered==0)~=0
    idx = filtered == 0;
    filtered(idx) = 0.0001;
end
rate = length(mask)./filtered;
rateDiff = imfilter(rate,grad);

if sum(isnan(rateDiff)) ~= 0
    idx = find(isnan(rateDiff));
    A = rateDiff(idx-1);
    B = rateDiff(idx +1);
    replacer = (B-A)/2+A;
    rateDiff(idx) = replacer;
end

dil = imdilate(rateDiff,SE);
startBurst = (dil - rateDiff) == 0;
%könnte noch durch einen noise filter ersetzt werden // statistik -->
%butterworth high pass, mean + 3*std = threshold;

%% filter parameters
interval = spikeTime(end)-spikeTime(1);
samples = length(rateDiff);
sampleRate = samples/interval;
lowPass = 0.024 * sampleRate;
highPass = 0.24 * sampleRate;
HPNoise = highPass/2;
settings.samplingFrequency = sampleRate;
settings.lowPass = lowPass;
settings.highPass = highPass;
settings.HPNoise = HPNoise;

%%
filterSettings  = constructFilt('Settings',settings);
[~,cleaner] =  filterSignal(rateDiff,filterSettings, 'high' );

starters = rateDiff >= 0.2*cleaner & startBurst;
%starters = rateDiff >=  startBurst;

rateDiffinv = -rateDiff;
dil = imdilate(rateDiffinv,SE);
stopBurst = (dil - rateDiffinv) == 0;
stops = rateDiffinv >= 0.2*cleaner & stopBurst;
%stops = rateDiffinv >= stopBurst;


%
idxStart = find(starters == 1);
idxStop = find(stops == 1)+1;

[idxStart, idxStop]=reSeq2(idxStart',idxStop',rate,0);


%% refinement
%[burstStart, burstStop] = refineEdges(idxStart, idxStop, spikeTime);
[burstStart, burstStop] = broadening_beta(idxStart, idxStop, spikeTime);

%[burstStart,burstStop] = spliceOver(burstStart,burstStop,length(spikeTime));
[burstStart, burstStop] = mergeBursts(burstStart, burstStop, spikeTime);

if burstStart(1) <= 0
    burstStart(1) = 1;
end




% merge if needed
%S = Maximum_median(abs(filtered)(1:end-length(mask)+2),40,'Type','percent','Dimension',2);
%IntraBurst = burstStart(1,2:end)-burstStop(1,1:end-1);

%% visualize
%visualizeBurst(rateDiff, burstStart, burstStop, spikeTime);



end



