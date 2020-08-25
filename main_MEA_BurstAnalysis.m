%UMGESTALTUNG --> BURSTING TOOL --> UNDER CONSTRUCTION

% This is a Burst Analysis Tool returning a bunch of novel burst parameter
% characterizing the phenotype of neurons cultured on MEA dishes
%----15.06.2020
%----written by Eva-Maria Weiss
%
% update: here data bins are read in instead the whole excel file list
%        containing the spike data
% main parts are therefore similar to 'MeaCalc'
%
% 
% (1) Detection of valid electrodes
% - for data cleaning
% - read in baseline data; this is basic for the whole experiment (also for
%   following days, each measurement that is linked to one experiment is
%   treated according to this result
%
% (2) Inititaliziation
%     data bins are loaded, cleared for invalid electrodes and stored as a
%     data cell per bin with timepoints x electrodes
%     a data struct is generated containing name of bin, data and the list
%     of valid electrodes
%
% (3) Grouping of Data per Well
    

function resultCell = main_MEA_BurstAnalysis(varargin)

%% set flags
%default

if nargin ~= 0    
    z = varargin{1};
end


%% (0) get all files
% be sure to have the files names in an ascending direction
[baseInfo, basePath] = uigetfile('*.csv','select data',...
                                 'MultiSelect','on');
[saveFile,savePath]= uiputfile('*.mat','save data');



%% (1) Data cleaning
%receive spike list of first baseline bin
if ~iscell(baseInfo)
    baseInfo = cellstr(baseInfo);
end

nbBin = length(baseInfo);

if ~exist('z','var')
    if length(baseInfo) >= 3
        z = 3;
    else
        z = 1;
    end
end

firstList = getList('instruction','insert baseline files','multiSelection','off',...
                    'folder',basePath,'file',baseInfo{z}); %receive the spike list of baseline
dataRaw = firstList.data;
[~,validList, validWells] = cleanData(dataRaw);


%% (2) Initializiation

dataStruct = repmat(struct('experiment',[]),1,nbBin);
dataStruct(1).validElectrods = validList;  
dataStruct(1).validWells = validWells;     

for bin = 1:nbBin
    %filename = strrep(baseInfo{bin},'.csv','');
    spikeList = getList('nameList',validList,'folder',basePath,'file',baseInfo{bin});
    % get time interval
    spikeListArray = spikeList.data;
    spikeVec = conversion(spikeListArray, 2); %get list of all time stamps from each electrode together
    interval = getInterval(spikeVec); %in seconds
    
    
    dataStruct(bin).experiment = spikeList.experiment{1,1};
    dataStruct(bin).data = spikeList.data;
    dataStruct(bin).interval = interval;
end

%% global parameter
wellNb = size(validWells,2);
resultCell = init_featureCell(wellNb,validWells,nbBin);
nbFeatures = size(resultCell,1)-1;

%%
for bin = 1:nbBin    
  
    spikeList = dataStruct(bin).data;
    interval = dataStruct(bin).interval; 
    [wellList, ~] = getWells(spikeList, -1, validWells);
    currNbWell = size(wellList,2);
    for well = 1:currNbWell    
        %check well
       
        numbSpikes = getSpikeNum(wellList{2,well});
        spikeRate = numbSpikes./interval;
        validnb = sum(spikeRate > 0.1);
        
        name = wellList(1,well);
        
        if validnb < 10 %miss number of valid electrodes condition
           wellIdx = find(strcmp(validWells,name{1,1}));
           resultCell(2:end,wellIdx+1,bin) = num2cell(nan);
        else        
        
        spikeTimes = wellList{2,well};
%sprintf('%i_%i',well,bin)
        burstVar = getNB_movFR(spikeTimes, interval);
        resultCell_onewell = featureCalculation(burstVar,interval,nbFeatures);
        wellIdx = find(strcmp(validWells,name{1,1}));
        resultCell(2:end,wellIdx+1,bin) = resultCell_onewell; 
        end
    end   
   
end
 filename = strrep(saveFile,'.mat','');
 save(fullfile(savePath, filename),'resultCell');

%% binning

end

