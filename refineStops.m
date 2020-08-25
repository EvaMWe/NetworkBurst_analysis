% this function will refine the start of a cluster by slope
% not used yet, does not seem to be necessary...
function StartRef = refineStops(stops,trace)
head = 100;
traceLong = [trace zeros(1,head)];
win = 2;
proc = 1;
%traceMat = repmat(traceLong,length(start),1);
adding = ones(1,length(stops));
nbSp = zeros(1,length(stops));
count = 0;
stop_temp = stops;
while proc ~= 0
    stop_temp = stop_temp + count;
    if traceLong(stop_temp(end)) == 0
        adding(end) = 0;
    end
    stopstop_temp = stop_temp + win;
    
    %section = traceMat(:,start+count-1:start+count-1+win);
    slope = arrayfun(@(a,b)  getSlope(a,b,traceLong), stop_temp, stopstop_temp);
    temp = slope < -10; %degree
    
    % no taking over
    logic = stop_temp(1,1:end-1) + count >= stop_temp(1,2:end);
    temp(logical([logic 1])) = 0;
%     
       
    adding = temp & adding; % save the zero!
    proc = sum(adding);
    nbSp = nbSp+adding;
    
    count = count+1;
end
StartRef = stops+nbSp;
end

% 
 function slope = getSlope(start,stop,trace)
 section = trace(start:stop);
 slope = calcSlopeNorm(section,trace, 'spline', 'degree');
 end