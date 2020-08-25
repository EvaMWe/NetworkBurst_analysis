%this function selects long bursts and clusters according todensity
%criterium 
% (1) calculate the density within detected bursts 
% (2) calculate the density outside the detected bursts

function [starts, stops] = selectBydense(start_, stop_, spiketimes,  varargin)

thresh = 0.1;

if nargin == 4
    thresh = varargin{1};
end

%% spike density within bursts
duration = spiketimes(stop_) - spiketimes(start_);
spNb = (stop_-start_)';
densIn = spNb./duration;

%% spike density outside bursts
spikeOut_1 =spiketimes(1,1:start_(1));
len1 = length(spikeOut_1);
spikeOut_2 =spiketimes(1,stop_(end):end);
len2 = length(spikeOut_2);

spikeOut = arrayfun(@(a,b) spiketimes(1,b:a),start_(2:end),stop_(1:end-1), 'UniformOutput',false);
len = length(spikeOut);
number = cellfun(@nnz, spikeOut);


spikeTimes_out = zeros(1,len1 + len2 + sum(number));

spNb = len1 + len2 + sum(number);
window = ceil(spNb/50);

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
S = (mean(densIn)-mean(densAll))/mean(densAll);

densMax = mean(densAll) + thresh*(max(densAll)-mean(densAll));


% 

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