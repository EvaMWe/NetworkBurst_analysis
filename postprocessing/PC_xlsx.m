% THis function will convert the data for PC1, calculating the mean,
% transform to a 2D array and export the data to excel

function PC_xlsx(resultCell)

%convert data to mean
nbBlocks = length(resultCell);

for bl = 1:nbBlocks
    dataBlock =resultCell{bl,1};
    nbGroups = size(dataBlock,3);
    nbTimepoints = size(dataBlock,2);
    meanCell = cell(nbTimepoints,nbGroups);
    
    for gr=1:nbGroups
        C = dataBlock(:,:,gr);
        M = cellfun(@(a) replace_empty(a), C);
        meanVec = mean(M,1,'omitnan');
        meanCell(1:length(meanVec),gr) = num2cell(meanVec);
    end
    
end
end

function r = replace_empty(a)
if isempty(a)
    r = nan;
else
    r = a;
end
end