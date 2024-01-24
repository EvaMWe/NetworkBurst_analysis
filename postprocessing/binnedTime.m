%this is a function to group the bin into experiment days

function resultCell_binned = binnedTime(varargin)

if nargin == 0
    [resultFile,resultPath] = uigetfile();
    resultCell = load(fullfile(resultPath,resultFile));
    resultCell=resultCell.resultCell;
else    
    resultCell = varargin{1};    
end
%% replace 0x0 double
wellNames = resultCell(1,2:end,1);
featureNames = resultCell(2:end,1);
resultCell_ = cellfun(@(x) checkEmpty(x),resultCell,'UniformOutput',false);

listOfGroups = xlsx2mat_grouping('sheetname', 'binningTime', 'nbGroups', 2,...
        'variableNames',{'Bin_start','Bin_stop'});
    start = cell2mat(listOfGroups(2:end,1));
    stop = cell2mat(listOfGroups(2:end,2));
    
    numResult = cell2mat(resultCell_(2:end,2:end,:));
    resultCell_binned = arrayfun(@(a,b) mean( numResult(:,:,a:b),3,'omitnan'), start, stop, 'UniformOutput', false);
    
    resultCell_binned{1,2} = wellNames;
    resultCell_binned{2,2} = featureNames;
end

function a = checkEmpty(a)
     if isempty(a) 
         a = 0;
     end
end