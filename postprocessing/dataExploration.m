% This is a function that will isolate the data for the relevant parameters
% and bring them to an exportable format
%Output file for each block of parameters containing the data for relevant
%features;
%for on parameters: time in row, well in columne,grouped (first wells of
%grouop1 than wells pf group2 ...)

% INPUT: result_fragments: n x m x h cell array containing data with n timepoints, m feature blocks; data in h=1; h=2 contains information about
% variable names (first row, for each block), well names(result_fragments_{1,2,2})

function resultCell = dataExploration(result_fragments, relevantFeatures, groupingCell)

[fileName, pathName] = uiputfile('*.xlsx','Select saving folder');


nbBlocks = size(result_fragments,2);
nbGroups = size(groupingCell,2);
nbTime = size(result_fragments,1);
nbWells = sum(sum(cellfun(@isscalar, groupingCell)));

for bl = 1:nbBlocks
    %get relevant features
    features = relevantFeatures{bl,1};
    nbF = cell2mat(features(:,2));
    nbFCum = cumsum(nbF);
    nbFrel = nbFCum./nbFCum(end);
    idx = find(nbFrel>0.8);
    nbF = idx(1);
    
    features = features(1:nbF,:);
       
    %get data block
    data = result_fragments(:,bl,1);
    varNames = result_fragments{bl,1,2};
    featureIndex = cellfun(@(nam) find(strcmp(varNames,nam)), features(:,1),'UniformOutput',false);
    wellNames = result_fragments{1,2,2};
    
    
    for f = 1:nbF
        fileName = extractBefore(fileName, '.');
        add = sprintf('FeatureBl_%i',bl);
        add = compose(add);
        fileName = strjoin(add,fileName);
        saveName = strcat(fullfile(pathName,fileName),'.xlsx');
        
        curFeature = featureIndex{f};
        sheetName = varNames{curFeature};
        if length(sheetName) > 31
            sheetName = sheetName(1:31);
        end
        resultCell = cell(nbTime+2,nbWells);
        count = 1;
        grIdx = cell(nbGroups,1);
        for group = 1:nbGroups
            wellNumber = groupingCell(:,group);
            idx = cellfun(@isscalar, wellNumber);             
            wellNumber = wellNumber(idx);
            grIdx{group,1} = wellNumber;
            
            resultCell(2,count:count+length(wellNumber)-1) = wellNames(cell2mat(wellNumber));
            resultCell{1,count} = groupingCell{1,group};
            count = count+length(wellNumber);
        end
        for t=1:nbTime
            data_t = data{t,1,1};
            data_f = data_t(curFeature,:);
            count = 1;
            for group = 1:nbGroups
                idx = cell2mat(grIdx{group});
                dataGr = data_f(idx);
                resultCell(t+2,count:count+length(dataGr)-1) = num2cell(dataGr);
                count = count+length(dataGr);
            end  
        warning('off','MATLAB:xlswrite:AddSheet');  
        sprintf('time:%i_block:%i_feature:%i',t,bl,f)
        writetable(cell2table(resultCell),saveName,'WriteVariableNames',false,...
        'Sheet',sheetName);
        end
        %export as excel file
       
    end
    
     % create result cell array
    
    
end

end
