% spikeTimes = cell array containing the spike times per electrode
% This is the main function for network burst analysiation according to the
% moving firinge rate method 
% Input (spikeTImes) is a cell array containing the data from one time bin
% derived from one network (column = electrodes, rows: spikeTimes (in
% seconds)
%
% This function is the core of burst detection
% (1)
% (2)
% (3) from 22.12.2020: combined detection of bursts getBurstCluster_IBI
% (based on smoothed IBI, using 1° derivative) and getBurstCluster_slope
% (based on the rate [ length(smoothingWindow) : IBI_smoothed]



function burstVar = getNB_movFR(spikeTimes, interval)
%%convert to vector
minNbContr = floor(size(spikeTimes,2)*0.25);
minBurstSpike = 3;

[spikeVec,~] = conversion(spikeTimes, 2);
spikeVec = sort(spikeVec);



%% BURST BINS + BURSTS
% (1) burst bins + burstlets
[burstletStart_, burstletStop_, burstBin_start_,burstBin_stop_] = getBurst_movDens_NB(spikeVec);
% (2) select true burst bins
[burstBinInfo, ~, burstBin_start, burstBin_stop] = burstExploration(burstBin_start_,burstBin_stop_,spikeVec,spikeTimes, minBurstSpike,minNbContr);
%% BURSTLETS
% (1) select true burstlets by conditions
[burstletInfo, ~, burstletStart, burstletStop] = burstExploration(burstletStart_, burstletStop_,spikeVec,spikeTimes, minBurstSpike,3);

%% NETWORK BURSTS
%-- get network bursts by merging burstlets--;
% (1) select burstlets from burstlets_ with at least 2 spike on 2 electrodes;
[~, ~, burstletcand_start, burstletcand_stop] = burstExploration(burstletStart_, burstletStop_,spikeVec,spikeTimes, 1,2);

% (2) form network bursts from burstlets
[burstStart_IBI,burstStop_IBI]= getBurstCluster_IBI(burstletcand_start, burstletcand_stop, spikeVec,0.3,0.3);
[burstStart_slope,burstStop_slope]= getBurstCluster_slope(burstletcand_start, burstletcand_stop, spikeVec,0.3,0.3);


[burstStart_, burstStop_] = combineBursts(burstStart_IBI,burstStop_IBI, burstStart_slope,burstStop_slope); 

[burstInfo, ~, burstStart_,burstStop_] = burstExploration(burstStart_,burstStop_,spikeVec,spikeTimes,minBurstSpike, minNbContr);
[burstStart,burstStop] = selectBydense_dyn(burstStart_,burstStop_, spikeVec,0.5);
[burstiStart_, burstiStop_] = getOutsideBursts(burstletcand_start, burstletcand_stop, burstStart,burstStop);

%burstlets outside of bursts
[burstOutInfo, ~, burstletOut_start, burstletOut_stop] = burstExploration(burstiStart_, burstiStop_,spikeVec,spikeTimes, minBurstSpike,3);


%%
% % CLUSTER
% % (5) burst clusters and burstlets outside clusters
% [Startcluster_, Stopcluster_, remainersStart,remainersStop] = findAgglomerats(burstletStart, burstletStop, lburstStart_,lburstStop_, spikeVec);
% if ~isempty(Startcluster_) && ~isempty(Stopcluster_)
%     [Startcluster, Stopcluster] = selectBydense(Startcluster_, Stopcluster_, spikeVec, 2);
%     [burstInfo_cl, ~, ~, ~] = burstExploration(Startcluster, Stopcluster,spikeVec,spikeTimes,minNbSpikes*2);
% else
%     Startcluster = [];
%     Stopcluster = [];
%     burstInfo_cl.nbcontribChannels = 0;
% end

% CLUSTER_PLUS_LONG_BURST_(CLB)
% (6) LongBurst outside cluster and cluster CLB
% [StartComplete_, StopComplete_] = selectLong(remainersStart, lburstStart, lburstStop);
% 
% if isempty(StartComplete_) || isempty(StopComplete_)
%     StartComplete =Startcluster;
%     StopComplete  = Stopcluster;
% elseif isempty(Startcluster) || isempty(Stopcluster)
%     StartComplete =StartComplete_;
%     StopComplete  = StopComplete_;
% else
%     StartComplete =sort([StartComplete_; Startcluster]);
%     StopComplete = sort([StopComplete_; Stopcluster]);
% end
% 
% if length(StartComplete) ~= length(StopComplete)
%     [StartComplete,StopComplete] = reSeq(StartComplete,StopComplete);
% end
% 
% [burstInfo_comp, ~, ~, ~] = burstExploration(StartComplete, StopComplete,spikeVec,spikeTimes,minNbSpikes*2);

%%


burstVar.spikeList = spikeVec;
burstVar.spikeData = spikeTimes;
burstVar.Interval = interval;
burstVar.burstBin_start = burstBin_start;
burstVar.burstBin_stop = burstBin_stop;
burstVar.burstletStart = burstletStart;
burstVar.burstletStop = burstletStop;
burstVar.NBstart= burstStart;
burstVar.NBstop = burstStop;
burstVar.burstletOut_start = burstletOut_start;
burstVar.burstletOut_stop= burstletOut_stop;

[burstVar.MContrCh_burstlet, burstVar.CV_ContrCh_burstlet, burstVar.Sk_ContrCh_burstlet] = statPara(extractfield(burstletInfo,'nbcontribChannels'));
[burstVar.MContrCh_burstBin, burstVar.CV_ContrCh_burstBin, burstVar.Sk_ContrCh_burstBin] = statPara(extractfield(burstBinInfo,'nbcontribChannels'));
[burstVar.MContrCh_NBburst, burstVar.CV_ContrCh_NBburst, burstVar.Sk_ContrCh_NBburst] = statPara(extractfield(burstInfo,'nbcontribChannels'));
[burstVar.MContrCh_NBburstOut, burstVar.CV_ContrCh_NBburstOut, burstVar.Sk_ContrCh_NBburstOut] = statPara(extractfield(burstOutInfo,'nbcontribChannels'));


end

