%this function selects long bursts and clusters according todensity
%criterium 
% (1) calculate the density within detected bursts 
% (2) calculate the density outside the detected bursts

function [starts, stops] = selectBydense(start_, stop_, spiketimes,  varargin)

thresh = 1.645; %is the 0.95 confidence interval
[start_,stop_] = reSeq(start_,stop_);
booster = 1.5;

if length(start_) <= 2 || length(stop_) < 2
    starts = start_;
    stops = stop_;
    return
end

if nargin == 4
   booster = varargin{1};
end

%% spike density within bursts
duration = spiketimes(stop_) - spiketimes(start_);
if size(stop_,1) > 1 
    stop_ = stop_';
end

if size(start_,1) > 1 
    start_ = start_';
end

spNb = (stop_-start_);
densIn = spNb./duration;
spNb_In = sum(spNb);

%% spike density outside bursts

start_ = start_';
stop_ = stop_';
spikeOut_1 =spiketimes(1,1:start_(1));
len1 = length(spikeOut_1);
spikeOut_2 =spiketimes(1,stop_(end):end);
len2 = length(spikeOut_2);


spikeOut = arrayfun(@(a,b) spiketimes(1,b:a),start_(2:end),stop_(1:end-1), 'UniformOutput',false);
len = length(spikeOut);
number = cellfun(@nnz, spikeOut);


spikeTimes_out = zeros(1,len1 + len2 + sum(number));

spNb_out = len1 + len2 + sum(number);
%% condition that at least a certain amount of spikes are outside bursts; otherwise threshold calculation will be invalide
if spNb_out/(spNb_out+ spNb_In) < 0.3
    starts = start_;
    stops = stop_;
    return
end

%%    
window = ceil(spNb_out/50);

spikeTimes_out(1,1:len1) = spikeOut_1;
count = len1;
for out = 1:len
    spikeTimes_out(1,count+1:count+number(out)) = cell2mat(spikeOut(out));
    count = count+number(out);
end
spikeTimes_out(1,end-len2+1:end) = spikeOut_2;

%split trace to calcualte threshold
linIndex = 1:window:length(spikeTimes_out);

durSplit = spikeTimes_out(linIndex(2:end)) -spikeTimes_out(linIndex(1:end-1));
densAll = window./durSplit;
if length(densAll) <= 10
    starts = start_;
    stops = stop_;
    return
end

%%
%S = (mean(densIn)-mean(densAll))/mean(densAll);

densMax = mean(densAll) + booster*thresh*(std(densAll)/(sqrt(length(densAll))));

idx = densIn > densMax;
%densAll =  [densIn densOut];
% densSorted = sort(densAll);
% densCum = cumsum(densSorted);
% densRel = densCum./densCum(end);
% 
% S = getSchnittpunkt(densRel);
% 
% densThresh = find(densRel > S);
% 
% 
% h = densSorted(densThresh(1));
% idx = densIn > h*theta;

starts = start_(idx);
stops = stop_(idx);

end