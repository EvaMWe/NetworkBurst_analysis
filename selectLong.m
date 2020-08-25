% this function is supposed to select long bursts that are concluded in
% remainers (bursts that are not inside clusters
%or better which remainer is a long burst

function [longStart, longStop] = selectLong(remainersStart, lburstStart, lburstStop)
count = 1;

if length(lburstStart) <= 1 || length(lburstStop) <= 1
    longStart = lburstStart;
    longStop = lburstStop;
    return
end

if lburstStop(1) < lburstStart(1)
    lburstStop(1) = [];
end

if length(lburstStart) ~= length(lburstStop)
    [lburstStart,lburstStop] = reSeq(lburstStart,lburstStop);
end

longStart = zeros(size(lburstStart));
longStop = zeros(size(lburstStart));

for k = 1:length(remainersStart)
    check = sum(lburstStart == remainersStart(k));
    if check ~= 0
        idx = find(lburstStart == remainersStart(k));
        longStart(count) = lburstStart(idx);
        longStop(count) = lburstStop(idx);
        count = count+1;
    end
end
longStart(longStart==0) = [];
longStop(longStop == 0) = [];

if length(longStart) ~= length(longStop)
    [longStart,longStop] = reSeq(longStart,longStop);
end

end