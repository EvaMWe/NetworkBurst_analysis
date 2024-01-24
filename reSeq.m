% this function iteratively takes start location and the next following
% stop location and stored in new vectors: finStop, finStart;
% Two following starts or stops are discarded in that way.

%input: (start,stop) 1xn double array containing all start and stop indices
%requirement: vectors are sorted from earliest to latest burst
%output: (finStart, finStop) 1xm double array with removed "extra" stops
%and starts

% PSEUDO
% finstart(1) = min(start) (= first value)
% finstop(1) = min(stop)

% count = 1
% WHILE no endpoint reached
%   finstart(count) = min(of all starts higher than finStop(count-1))
%   finstop(count) = min(of all stops higher than finStart(count-1))
% END

function [finStart,finStop]=reSeq(start,stop)

if isempty(start) || isempty(stop)
    finStart = [];
    finStop = [];
    return
end

if size(start,2) > 1
   start = start';
end

if size(stop,2) >1
    stop = stop';
end

if length(start) < 2     
    if start(1) < stop(1) 
    finStart = start(1);
    finStop = stop(1);
    elseif start > stop(1) && length(stop) > 1
        finStart = start(1);
        finStop = stop(2);
    else
        finStart = [];
        finStop = [];
    end
        
    return
end

if  length(stop) < 2
    finStop = stop;
    finStart = start(1);
    if start(1) > stop(1)
        finStart = [];
        finStop = [];
    end  
    
    return
end

if start(1) > stop(1)
    stop=stop(2:end,1);
end



endStart = stop(end,1);
endStop = start(end,1);

finStart = zeros(length(start),1);
finStop = zeros(length(stop),1);

%initialize
finStart(1,1) = min(start);
finStop(1,1) = min(stop);


count = 2;
cont = 1;
while cont == 1
%sprintf('%i',count)
   if length(finStop) < count || length(finStart) < count
       cont = 0;
        continue
   elseif finStop(count) >= endStop || finStart(count) >= endStart || finStart(count) >= endStop || finStop(count) >= endStart
        cont = 0;
        continue
   end
    
    tail_start = start(start>finStop(count-1));        
    if ~isempty(tail_start) 
        finStart(count) = min(tail_start);        
    else
        cont = 0;
        continue
    end
    
    tail_stop = stop(stop>finStart(count));
    if  ~isempty(tail_stop)        
        finStop(count) = min(tail_stop);
    else
        cont = 0;
        continue
    end
 
    count = count+1;
end

finStart(finStart == 0) =[];
finStop(finStop == 0)=[];

%check for equal length
if length(finStart) > length(finStop)
    finStart = finStart(1:length(finStop));
end

if length(finStop) > length(finStart)
    finStop = finStop(1:length(finStart));
end

% check if start value is lower than first stop value // could be changed
% due to initial problem (start(1,1) gets start_temp(1,1) if start(1,1) is
% in the negative range
if finStart(1,1) > finStop(1,1)
    finStart = finStart(2:end,1);
    finStop = finStop(2:end,1);
end


end
