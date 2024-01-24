%this function selects long bursts and clusters according todensity
%criterium 
% (1) calculate the density within detected bursts 
% (2) calculate the density outside the detected bursts

function [starts, stops] = selectBydense_dyn(start_, stop_, spiketimes,  varargin)


[start_,stop_] = reSeq(start_,stop_);

if length(start_) <= 2 || length(stop_) < 2
    starts = start_;
    stops = stop_;
    return
end

booster = 1.5;
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

%% spike density outside bursts

start_ = start_';
stop_ = stop_';
spikeOut_1 =spiketimes(1,1:start_(1));
len1 = length(spikeOut_1);
spikeOut_2 =spiketimes(1,stop_(end):end);
len2 = length(spikeOut_2);


spikeOut = arrayfun(@(a,b) spiketimes(1,b:a),start_(2:end),stop_(1:end-1), 'UniformOutput',false);

number = cellfun(@nnz, spikeOut);
interval = cellfun(@(a) a(end)-a(1), spikeOut);
dense = number./interval;
dense_1 = len1/(spikeOut_1(end)-spikeOut_1(1));
dense_2 = len2/(spikeOut_2(end)-spikeOut_2(1));
dense = [dense_1 ; dense ; dense_2];
dense_smth = smooth(dense,3);
thresh = (dense_smth.*booster)';



idx = densIn > thresh(1,1:length(densIn));


starts = start_(idx);
stops = stop_(idx);

end