% this function will refine the start of a cluster by slope
% not used yet, does not seem to be necessary...
function StartRef = refineStarts(start,trace)
head = 100;
traceLong = [zeros(1,head) trace];
win = 2;
proc = 1;
%traceMat = repmat(traceLong,length(start),1);
adding = ones(length(start),1)';
nbSp = zeros(1,length(start));
count = 0;
stop_temp = start+head;
while proc ~= 0
    stop_temp = stop_temp - count;
    if traceLong(stop_temp(1)) == 0
        adding(1) = 0;
    end
    start_temp = stop_temp - win;
    %section = traceMat(:,start+count-1:start+count-1+win);
    slope = arrayfun(@(a,b)  getSlope(a,b,traceLong), start_temp, stop_temp);
    temp = slope > 20; %degree
    
    % no taking over
    logic = start_temp(1,2:end) - count <= start_temp(1,1:end-1);
    temp(logical([1 logic])) = 0;
%     
       
    adding = temp & adding; % save the zero!
    proc = sum(adding);
    nbSp = nbSp+adding;
    
    count = count+1;
end
StartRef = start-nbSp;
end

% 
 function slope = getSlope(start,stop,trace)
 section = trace(start:stop);
 slope = calcSlopeNorm(section,trace, 'spline', 'degree');
 end