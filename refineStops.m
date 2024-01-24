% this function will refine the start of a cluster by slope
% not used yet, does not seem to be necessary...
% 0 = no dynamic threshold
% 1 = dynamic threshold
% Note: dynamic threshold means a vector of the same size as trace
% containing corresponding threshold values for each point

function StartRef = refineStops(stops,trace,h, varargin)
head = 100;
dyn = 0;

if nargin == 4
    dyn = varargin{1};
    h_ = h;
    h_ = [h_ repmat(h_(end),1,head)];
end

traceLong = [trace repmat(-1,1,head)];
win = 2;
proc = 1;
%traceMat = repmat(traceLong,length(start),1);
adding = ones(1,length(stops));
nbSp = zeros(1,length(stops));
count = 0;
stop_temp = stops;
while proc ~= 0    
    if traceLong(stop_temp(end)) == -1
        adding(end) = 0;
    end
    stopstop_temp = stop_temp + win;
    
    %section = traceMat(:,start+count-1:start+count-1+win);
    slope = arrayfun(@(a,b)  getSlope(a,b,traceLong), stop_temp, stopstop_temp);
    if dyn == 1
        h = h_(stop_temp);
    end
    temp = slope < -10 | traceLong(stop_temp) > h; %degree
    
    % no taking over
    logic = stop_temp(1,1:end-1) + count >= stop_temp(1,2:end);
    temp(logical([logic 0])) = 0;
%     
       
    adding = temp & adding; % save the zero!
    proc = sum(adding);
    nbSp = nbSp+adding;
    stop_temp = stop_temp + 1;
    count = count + 1;
end
StartRef = stops+nbSp;
end

% 
 function slope = getSlope(start,stop,trace)
 section = trace(start:stop);
 slope = calcSlopeNorm(section,trace, 'spline', 'degree');
 end