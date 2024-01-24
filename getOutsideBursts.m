function [burstOutStart, burstOutStop] = getOutsideBursts(allB_starts, allB_stops, LB_start, LB_stop)
 Bstart = repmat(allB_starts,length(LB_start),1);
    Log1 = Bstart < LB_start;
    Log2 = Bstart > LB_stop;
    Log = Log1 + Log2;
    sumLog = sum(Log,1);
    shStartlog = sumLog == length(LB_start);
    burstOutStart = allB_starts(shStartlog);
    % get associated stop
    burstOutStop = allB_stops(shStartlog);
end