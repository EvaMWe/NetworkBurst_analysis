% this is a function to average bins from the same experimental phase
% resultCell contains the features per well and bin, while bins are in
% direction 3.
% timepoints contains the start and stop bin number; it's the frame within
% the bin are summarized by calculating the mean; it is a double array;
% first column: starts; second column: stops;

function resultCell_binned = binning(resultCell,timepoints)
start = timepoints(:,1);
stop = timepoints(:,2);
resultCell_array = resultCell(2:end,2:end,:);
resultCell_binarray = zeros(size(resultCell_array,1),size(resultCell_array,2), length(timepoints))

resultCell_binned = arrayfun(@(a,b) mean(resultCell(:,:,a:b),3),start,stop);
end