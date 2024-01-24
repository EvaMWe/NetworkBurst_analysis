% this function take the task of reSeq:
% this function is used in concert with the getBurst_movDens_NB function

%start,stop = 1xn double
%gradStart, gradStop = 1xn double

% this function will give the stop indices a higher priority as being
% fixed since they have been validated befor

function [Start, Stop]= reSeq3(start,stop,trace, varargin)
%default
recalc = 0;

if nargin == 5
    recalc = varargin{1};
end

if length(start) <=2 || length(stop) <=2
    Start = start;
    Stop = stop;
    return
end

if size(start,1) > 1
    start = start';
end

if size(stop,1) > 1
    stop = stop';
end

while stop(1) < start(1)
    stop = stop(2:end);
end
lenstop = length(stop);


% falls stops gesucht werden


P =repmat(start',1,lenstop);

logi = P < stop;

difference = diff(sum(logi));

%% get indices of extra starts

checkStarts = diff(difference == 2);
if sum(checkStarts) == 0
    Start = start;
else
    
    e1_c = find(checkStarts == 1)+1;
    e2_c = find(checkStarts == -1)+1;
    
    e1 = sum(logi(:,e1_c));
    e2 = sum(logi(:,e2_c));
    if length(e1) < 2 || length(e2) < 2
        Start = start;
    else
        
        if e1(1) >= e2(1)
            e2 = e2(2:end);
        end
        
        if e1(end) >= e2(end)
            e1 = e1(1:end-1);
        end
        
        startFin = e1;
        prop_col_start = difference < 2;
        properStartsidx = sum(logi(:,prop_col_start));
        properStart = start(properStartsidx);
        % add last
        if difference(end) == 1
        lasti = start(sum(logi(:,end)));
        properStart = [properStart lasti];
        end
        Start = unique(sort([properStart startFin]));
        %% get proper spikes
        
    end
end

checkStops = difference == 0;
if sum(checkStops) == 0
    Stop = stop;
else
    
    checkStops   = diff(checkStops);
    e1 = find(checkStops == 1)+1;
    e2 = find(checkStops == -1)+1;
    
    if length(e1) < 2 || length(e2) < 2
        Stop = stop;
    else
        
        if difference(end) == 0
            e2 = [e2 length(difference)+1];
        end
        
        if e1(1) > e2(1)
            e2 = e2(2:end);
        end
        
        if e1(end) > e2(end)
            e1 = e1(1:end-1);
        end
        
        
        stopFin = arrayfun(@(a,b) getMin(a,b,stop,trace,recalc), e1, e2);
        prop_col_stops = difference > 0;
        properStops = stop(prop_col_stops);
        %add last value:
        if difference(end) == 1
            properStops = [properStops stop(end)];
        end
        
        %% get together
        
        Stop = unique(sort([properStops stopFin]));
    end
end

if Start(end) > Stop(end)
    Start(end) = [];
end
Stop(Stop <= Start(1)) = [];
[Start,Stop]=reSeq(Start,Stop);


end

function Fin = getMin(idx1,idx2,spikes,trace, recalc)
segment = trace(spikes(idx1):spikes(idx2));
% recalculate idx 2 if a zero crossing is inbetween

if sum(segment <= 0) && recalc == 1
    spike_zero_ = find(segment<=0);
    spike_zero = spike_zero_(1)+spikes(idx1)-1;
    idx2_ = find(trace(spikes(idx1):spikes(idx2)) <= spike_zero);
    idx2 = idx2_(1)+idx1-1;
end

[~,Fin] = min(trace(spikes(idx1:idx2)));
Fin = spikes(idx1+Fin-1);
end
