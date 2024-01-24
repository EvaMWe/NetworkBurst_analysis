%this function helps to find the number of strings in a cell array of tpye
%string

function numberOfStrings = countStrings(string_ori,varargin)

if nargin == 1
    sorting = 0;
elseif nargin == 2
    sorting =varargin{1,1};
else
    error('invalid number of inputs')
end

stringNames = unique(string_ori);

numberOfStrings = cell(length(stringNames),2);
numberOfStrings(:,1) = stringNames;
stringnumber = arrayfun(@(nam) nnz(strcmp(string_ori,nam)), stringNames);
numberOfStrings(:,2) = num2cell(stringnumber);

if sorting == 1
    [~,sortedIdx]= sort(stringnumber,1,'descend');
    numberOfStrings_ = numberOfStrings;
    numberOfStrings = numberOfStrings_(sortedIdx,:);
end
    
    
end
