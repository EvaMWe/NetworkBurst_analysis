


function groupingCell= groupNBAnalysis_forNB(dataCell,varargin)

if nargin == 1
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-firstGrouping', 'nbGroups', 10 ,'variableNames',{'m1','m2','m3','m4','m5','m6','m7','m8','m9','m10'});
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-firstGrouping', 'nbGroups', 7 ,...
        %'variableNames',{'m1','m2','m3','m4','m5','m6','m7'});
   listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-groupedPerGene', 'nbGroups', 2 ,'variableNames',{'wT','SHP2'});
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-groupedPerGene', 'nbGroups', 2 );
elseif strcmp(varargin{2},'table parameters')
    sheetname = varargin{3};
    nbGroups = varargin{4};
    variableNames = varargin{5};
    listOfGroups = xlsx2mat_grouping('sheetname', sheetname, 'nbGroups', nbGroups ,'variableNames',variableNames);
else
    listOfGroups = varargin{2};
end


% get the cell array out of the data struct
wellNames = dataCell{1,2,2};
nbG = size(listOfGroups,2);
maxNbWell = size(listOfGroups,1)-1;

groupingCell = cell(maxNbWell+1,nbG);
groupingCell(1,:) = listOfGroups(1,:);

for gr = 1:nbG
    names = listOfGroups(2:end,gr);
    numbered = arrayfun(@(a) find(strcmp(wellNames,a)), names, 'UniformOutput',false);
    numbered = cellfun(@(a) replace(a), numbered);    
    numbered = sort(numbered);
    numbered(numbered == 0) = [];
    len = length(numbered);
    groupingCell(2:len+1,gr) = num2cell(numbered);
end

end

function r = replace(a)
if isempty(a)
    r = 0;
else
    r = a;
end
end
    
