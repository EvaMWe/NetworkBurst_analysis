

function stat_burstCalculator = groupNBfeatures(varargin)

if nargin == 0
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-firstGrouping', 'nbGroups', 4 ,'variableNames',{'m1','m2','m3','m4'});
    listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-firstGrouping', 'nbGroups', 10 ,'variableNames',{'m1','m2','m3','m4','m5','m6','m7','m8','m9','m10'});
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-groupedPerGene', 'nbGroups', 2 ,'variableNames',{'wT','SHP2'});
    %listOfGroups = xlsx2mat_grouping('sheetname', 'SHP2-groupedPerGene', 'nbGroups', 2 )
else
    listOfGroups = varargin{1};
end

[fileName, pathName] = uigetfile('Select files for grouping', 'MultiSelect','on');

% get the cell array out of the data struct
nbM = length(fileName);
nbG = size(listOfGroups,2);
maxNbWell = size(listOfGroups,1)-1;

if ~iscell(fileName)
    fileName = cellstr(fileName);
end

 %preallocate datastruct stat_spikeCalculator
 stat_burstCalculator = repmat(struct('measurement_name',[]),1,nbM);
for m = 1:nbM
    dataStruct = load(fullfile(pathName,fileName{m}),'-mat');
    data = dataStruct.burstCalcResult.summary_of_Results;
    namWell = data(1,2:end);
    nbWell = length(namWell);
    
    %preallocate cells for data storage // one cell per parameter
    numberOfBursts = cell(maxNbWell+2, nbG);
    MBR = cell(maxNbWell+2, nbG); %mean burst rate
    wMBR = cell(maxNbWell+2, nbG);
    MBD = cell(maxNbWell+2, nbG); %mean burst duration
    SDBD = cell(maxNbWell+2, nbG); %standard deviation of burst duration
    avg_nbSpikesInBurst = cell(maxNbWell+2, nbG); 
    std_nbSpikesInBurst = cell(maxNbWell+2, nbG);
    
   
    for gr = 1:nbG
        grNam = listOfGroups{1,gr};
        wells = listOfGroups(2:end,gr);     
        
        
        %nbVal = length(wells);
        dataSub = data(2:end,2:end);
        log = zeros(nbWell,1);
        for k = 1:nbWell
            log(k) = sum(strcmp(namWell{k},wells)) >= 1;
        end
        log = logical(log);
        
        %number of bursts per well       
        numberOfBursts(1,1) = {'number of bursts per well'};
        numberOfBursts(2,gr) = {grNam};
        numberOfBursts(3:length(dataSub(1,log')')+2,gr) = dataSub(1,log')';
        
        %mean burst rate per well     
        MBR(1,1) = {'mean burst rate per well'};
        MBR(2,gr) = {grNam};
        MBR(3:length(dataSub(1,log')')+2,gr) = dataSub(2,log')';
        
        
        %weighted mean burst rate per well     
        wMBR(1,1) = {'weighted mean burst rate per well'};
        wMBR(2,gr) = {grNam};
        wMBR(3:length(dataSub(1,log')')+2,gr) = dataSub(3,log')';
        
        %mean burst duration per well
        MBD(1,1) = {'mean burst duration per well'};
        MBD(2,gr) = {grNam};
        MBD(3:length(dataSub(1,log')')+2,gr) = dataSub(4,log')';
        
        %standard deviation of burst duration per well
        SDBD(1,1) = {'standard deviation of burst duration per well'};
        SDBD(2,gr) = {grNam};
        SDBD(3:length(dataSub(1,log')')+2,gr) = dataSub(5,log')';
        
        %average number of spikes per burst
        avg_nbSpikesInBurst(1,1) = {'average number of spikes per burst'};
        avg_nbSpikesInBurst(2,gr) = {grNam};
        avg_nbSpikesInBurst(3:length(dataSub(1,log')')+2,gr) = dataSub(6,log')';
        
        %standard deviation of number of spikes per burst
        std_nbSpikesInBurst(1,1) = {'standard deviation of number of spikes per burst'};
        std_nbSpikesInBurst(2,gr) = {grNam};
        std_nbSpikesInBurst(3:length(dataSub(1,log')')+2,gr) = dataSub(6,log')';
        
        go2stat_burstCalculator.numberOfBursts = numberOfBursts;
        go2stat_burstCalculator.numbContrEl = MBR;
        go2stat_burstCalculator.wMBR = wMBR;
        go2stat_burstCalculator.MBD = MBD;
        go2stat_burstCalculator.SDBD = SDBD;
        go2stat_burstCalculator.mean_nbSpikesBurst = avg_nbSpikesInBurst;
        go2stat_burstCalculator.std_nbSpikesBurst = std_nbSpikesInBurst;
        
    end
    stat_burstCalculator(m).measurement_name = strrep(fileName{m},'.mat','');
    stat_burstCalculator(m).data_for_statistic = go2stat_burstCalculator;
end
savename = fullfile(pathName,'BurstCalculator_Evaluation');
if ~exist (savename,'dir')
    mkdir(pathName,'BurstCalculator_Evaluation');
end
save (fullfile(savename,'burstCalc_grouped'), 'stat_burstCalculator');

end