% this function will refine the start of a cluster by slope
% not used yet, does not seem to be necessary...
% 0 = no dynamic threshold
% 1 = dynamic threshold
% Note: dynamic threshold means a vector of the same size as trace
% containing corresponding threshold values for each point
function StartRef = refineStarts(start,trace,h, varargin)
head = 100;
dyn = 0;
if nargin == 4
    dyn = varargin{1};
    h_ = h;
    h_ = [repmat(h_(1),1,head) h_];
end

traceLong = [zeros(1,head) trace];
win = 2;
proc = 1;
%traceMat = repmat(traceLong,length(start),1);
adding = ones(length(start),1)';
nbSp = zeros(1,length(start));
count = 0;
stop_temp = start+head;
while proc ~= 0
    
    if traceLong(stop_temp(1)) == 0
        adding(1) = 0;
    end
    start_temp = stop_temp - win;
    %section = traceMat(:,start+count-1:start+count-1+win);
    slope = arrayfun(@(a,b)  getSlope(a,b,traceLong), start_temp, stop_temp);
    if dyn == 1
        h = h_(stop_temp);
    end
    temp = slope > 20 | traceLong(stop_temp) > h; %degree
    
    % no taking over
    logic = start_temp(1,2:end) - count <= start_temp(1,1:end-1);
    temp(logical([0 logic])) = 0;
%     
       
    adding = temp & adding; % save the zero!
    
    proc = sum(adding);
    nbSp = nbSp+adding;
    stop_temp = stop_temp - 1;
    count = count+1;
end
StartRef = start-nbSp;
end

% 
 function slope = getSlope(start,stop,trace)
 section = trace(start:stop);
 slope = calcSlopeNorm(section,trace, 'spline', 'degree');
 end