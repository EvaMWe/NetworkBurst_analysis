function [resultTable_PC1, relevantFeatures, groupingCell] = PCA_Routine(result_fragments)
nbBlocks = size(result_fragments,2);
timepoints = size(result_fragments,1);

groupingCell= groupNBAnalysis_forNB(result_fragments);
nbGroups = size(groupingCell,2);

relevantFeatures = cell(nbBlocks,1);
wellNames = result_fragments{1,2,2};
resultTable_PC1 = cell(nbBlocks,1);
for bl = 1:nbBlocks
    varNames = result_fragments{bl,1,2};
    coursePC1 = cell(length(groupingCell), timepoints, nbGroups);
    relFeatures = cell(timepoints*length(varNames)^2,1);
    cnt = 1;
    for t= 1:timepoints
        data = result_fragments{t,bl,1};
        %% cleaning and converting
        % nan to 0
       
        data(isnan(data)) = 0;
        idxZero = find(any(data,1));
        validWells = wellNames(idxZero);
        data_ = data';
        [featureVector,PC, eigVal] = PCA_MEA_cell(data_);
        %% get most relevant features for that timepoint and block
        if isempty(PC) 
            continue
        end
        
        relFeatures_ = getRelevants(featureVector,eigVal,varNames);
        relFeatures(cnt:cnt+length(relFeatures_)-1) = relFeatures_;
        cnt = cnt+length(relFeatures);
        
        idxCell = ~(cellfun(@isempty, relFeatures));
        relFeatures = relFeatures(idxCell);

       
        
        %grouping 
        for gr = 1:nbGroups
            sprintf('%i_%i_%i',bl,t,gr)
            indexGr = cell2mat(groupingCell(2:end,gr));
            log = arrayfun(@(a) ismember(a,idxZero), indexGr);
            indexValid = indexGr(log); 
            PC1 = PC(indexValid,1);
            coursePC1(1:length(PC1),t,gr) = num2cell(PC1);
        end
        
        
        %fraction = proceedPC(PC,eigVal, feautureVector,groupingCell,idxZero);
    end
    numberOfRelevants = countStrings(relFeatures,1);
    relevantFeatures{bl,1} = numberOfRelevants;
    resultTable_PC1{bl,1} = coursePC1;
end
