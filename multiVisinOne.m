%Tool to visualize burst detection
%spikeTimes : cell array containing spike data per electrode

%---HINTS
%spikeSamples: spikeTimes converted to sample points by *12500 

%varargin{1}:
% type: 1 = input is time
%       2 = input is spike stamp (iinteger)

function multiVisinOne(spikeTimes,burstStart,burstStop, varargin)
type = 1;

if nargin == 4
    type = varargin{1};
end

spikeArray = conversion(spikeTimes,1);
nbEl =size(spikeArray,2);

spikeVec = sort(conversion(spikeTimes,2));

firstPos = min(spikeVec)*12500;
lastPos =max(spikeVec)*12500;



fig = figure ('NumberTitle', 'off', 'Name', 'timestamp plot');
place = 0.3;
for pl = 1:nbEl
    hold on
    spikeSamples = ceil(spikeArray(:,pl)*12500);
    spikeSamples(spikeSamples == 0) = [];
%     spikeStemps = zeros(nbSamples,1);
%     spikeStemps(spikeSamples) = 1;
%     
    %% create line coordinates
    spikeSamples_ = spikeSamples';
    x = repmat(spikeSamples_,2,1);
    y = repmat(pl-1,2,length(spikeSamples));
    y(2,:) = y(2,:) + 1 - place;
    
    line(x,y, 'color','k');

end

if type == 2
    burstStempStart = burstStart;
    burstStempStop = burstStop;
    
else
    burstStempStart = floor(spikeVec(burstStart).*12500);
    burstStempStop = ceil(spikeVec(burstStop).*12500);
end

width = (burstStempStop-burstStempStart);
widthVal = width(width > 0); 
x = burstStempStart(width>0);
nbBurst = length(x);

y = repmat(-0.2,1,nbBurst);
height = repmat(nbEl + 0.4,1,nbBurst);
for b = 1:nbBurst    
    rectangle('Position', [x(b),y(b),widthVal(b), height(b)], 'EdgeColor', 'red');
end

end

