function resultCell_ = featureCalculation(burstVar,interval,nbFeatures)

%%
resultCell_  = cell(nbFeatures,1);
spikeList = burstVar.spikeList;

nbSpikes = length(spikeList); %number of spikes
resultCell_{1,1} = nbSpikes;  
resultCell_{2,1} = nbSpikes/interval; %mean firing rate
    ISI_s = diff(spikeList);
[resultCell_{3,1}, resultCell_{4,1}, resultCell_{5,1}] = statPara(ISI_s);

burstletInvalid = 0;
burstInvalid = 0;


%% BURSTBINS
nbBins = length(burstVar.burstBin_start);
resultCell_{6,1} = nbBins;
resultCell_{7,1} = nbBins/interval;
resultCell_{8,1} = burstVar.MContrCh_burstBin;
resultCell_{9,1} = burstVar.CV_ContrCh_burstBin;
resultCell_{10,1} = burstVar.Sk_ContrCh_burstBin;


%% BURSTLETS
burstletStart = burstVar.burstletStart;
burstletStop = burstVar.burstletStop;
nbBurstlets = length(burstletStart);

if ~isempty(nbBurstlets) && ~isempty(burstletStop)&& ~isempty(burstletStart)
    resultCell_{11,1} = nbBurstlets;
    resultCell_{12,1} = nbBurstlets/interval;
    %IBI
    IBI_burstlets  = spikeList(burstletStart(2:end)) - spikeList(burstletStop(1:end-1));
    [resultCell_{13,1},resultCell_{14,1},resultCell_{15,1}] = statPara(IBI_burstlets);
    %duration
    D_burstlets =spikeList(burstletStop) - spikeList(burstletStart);
    [resultCell_{16,1},resultCell_{17,1},resultCell_{18,1}] = statPara(D_burstlets);
    %number of spikes in bursts
    nbS_burstlets = burstletStop - burstletStart; 
    [resultCell_{19,1},resultCell_{20,1},resultCell_{21,1}] = statPara(nbS_burstlets); 
    %firing rate
    FR_burstlets = nbS_burstlets/D_burstlets;
    [resultCell_{22,1},resultCell_{23,1},resultCell_{24,1}] = statPara(FR_burstlets);   
    %MFR outside burstlets
    nbS_outsideBurstlets =  burstletStart(2:end)-1 - burstletStop(1:end-1)+1;
    resultCell_{25,1} = (burstletStart(1)-1 + sum(nbS_outsideBurstlets) + length(spikeList)-burstletStop(end)+1)/interval; %add tails
    [outsideSpikes] = isolate_oSpikes(burstletStart,burstletStop, spikeList);
    ISI_outBurstlets =  diff(outsideSpikes);
    [resultCell_{26,1},resultCell_{27,1},resultCell_{28,1}] = statPara(ISI_outBurstlets);
    resultCell_{29,1} = burstVar.MContrCh_burstlet;
    resultCell_{30,1} = burstVar.CV_ContrCh_burstlet;
    resultCell_{31,1} = burstVar.Sk_ContrCh_burstlet;
    
else
    burstletInvalid = 1;
end

NBstart = burstVar.NBstart;
NBstop = burstVar.NBstop;
nb_bursts = length(NBstart);

if length(NBstart)>2 && length(NBstop)>2 && nb_bursts>2
    MBR_bursts = nb_bursts/interval;
    resultCell_{32,1} = nb_bursts;
    resultCell_{33,1} = MBR_bursts;
    %inter burst intervals
    IBI_bursts  = spikeList(NBstart(2:end)) - spikeList(NBstop(1:end-1));
    [resultCell_{34,1},resultCell_{35,1},resultCell_{36,1}] = statPara(IBI_bursts);
    %mean duration
    D_bursts =spikeList(NBstop) - spikeList(NBstart);
    [resultCell_{37,1},resultCell_{38,1},resultCell_{39,1}] = statPara(D_bursts);
    
    %number of spikes in bursts
    nbS_bursts= NBstop - NBstart; 
    [resultCell_{40,1},resultCell_{41,1},resultCell_{42,1}] = statPara(nbS_bursts); 
    %firing rate
    if burstletInvalid ~= 1
        FR_bursts = nbS_burstlets/D_burstlets;
        [resultCell_{43,1},resultCell_{44,1},resultCell_{45,1}] = statPara(FR_bursts);
        [outsideSpikes] = isolate_oSpikes(burstletStart,burstletStop, spikeList);
        ISI_outBursts =  diff(outsideSpikes);
        [resultCell_{47,1},resultCell_{48,1},resultCell_{49,1}] = statPara(ISI_outBursts);
    end
    %MFR outside bursts
    nbS_outsideBursts =  NBstart(2:end)-1 - NBstop(1:end-1)+1;
    resultCell_{46,1} = (NBstart(1)-1 + sum(nbS_outsideBursts) + length(spikeList)-NBstop(end)+1)/interval; %add tails
    %ISI outside bursts
   
else
    burstInvalid = 1;
end

if burstInvalid==0 && burstletInvalid==0
    [burstletinside_start, burstletinside_stop, nbburstlet_NB] = selectShinLong(burstletStart, burstletStop, NBstart, NBstop);
    [resultCell_{50,1},resultCell_{51,1},resultCell_{52,1}] = statPara(nbburstlet_NB);   
    %duration of burstlets inside NB
    D_burstlet_NB = spikeList(burstletinside_start) - spikeList( burstletinside_stop);
    [resultCell_{53,1},resultCell_{54,1},resultCell_{55,1}] = statPara(D_burstlet_NB);
    %interspike interval of burstlets inside NBs
    IBI_burstlets_in_NB  = spikeList(burstletinside_start(2:end)) - spikeList(burstletinside_stop(1:end-1));
    [resultCell_{56,1},resultCell_{57,1},resultCell_{58,1}] = statPara(IBI_burstlets_in_NB);
    %----outside----   
    burstletOut_start =burstVar.burstletOut_start;
    burstletOut_stop =burstVar.burstletOut_stop;
    % mean bursting rate of burstlets outside NB
    resultCell_{59,1} = length( burstletOut_start)/(interval-sum(D_burstlet_NB));
    % interbursting interval of  burstlets outside NB
    IBI_burstlets_outNB  = spikeList(burstletOut_start(2:end)) - spikeList(burstletOut_stop(1:end-1));
    [resultCell_{60,1},resultCell_{61,1},resultCell_{62,1}] = statPara(IBI_burstlets_outNB);
    resultCell_{63,1} = burstVar.MContrCh_NBburst;
    resultCell_{64,1} = burstVar.CV_ContrCh_NBburst;
    resultCell_{65,1} = burstVar.Sk_ContrCh_NBburst;
    resultCell_{66,1} = burstVar.MContrCh_NBburstOut;
    resultCell_{67,1} = burstVar.CV_ContrCh_NBburstOut;
    resultCell_{68,1} = burstVar.Sk_ContrCh_NBburstOut;    
end
%%

end