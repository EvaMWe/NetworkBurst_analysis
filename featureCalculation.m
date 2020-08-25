function resultCell_ = featureCalculation(burstVar,interval,nbFeatures)

%%
resultCell_  = cell(nbFeatures,1);
spikeList = burstVar.spikeList;

nbSpikes = length(spikeList);
resultCell_{1,1} = nbSpikes;
resultCell_{2,1} = nbSpikes/interval;
    ISI_s = diff(spikeList);
[resultCell_{3,1}, resultCell_{4,1}, resultCell_{5,1}] = statPara(ISI_s);

shortInvalid = 0;
longInvalid = 0;
clusterInvalid = 0;
compInvalid = 0;

%% 
SB_start = burstVar.shortBurst_starts;
SB_stop = burstVar.shortBurst_stops;
nbSB = length(SB_start);
if ~isempty(nbSB) && ~isempty(SB_stop)&& ~isempty(SB_start)
    resultCell_{6,1} = nbSB;
    resultCell_{7,1} = nbSB/interval;
    ISBI_s  = spikeList(SB_start(2:end)) - spikeList(SB_stop(1:end-1));
    [resultCell_{8,1},resultCell_{9,1},resultCell_{10,1}] = statPara(ISBI_s);
    DSB =spikeList(SB_stop) - spikeList(SB_start);
    [resultCell_{11,1},resultCell_{12,1},resultCell_{13,1}] = statPara(DSB);
    nbS_SB = SB_stop - SB_start;
    [resultCell_{14,1},resultCell_{15,1},resultCell_{16,1}] = statPara(nbS_SB);
    nbS_oSB =  SB_start(2:end)-1 - SB_stop(1:end-1)+1;
    resultCell_{17,1} = SB_start(1)-1 + sum(nbS_oSB) + length(spikeList)-SB_stop(end)+1;
    [outsideSpikes] = isolate_oSpikes(SB_start,SB_stop, spikeList);
    ISI_oSB =  diff(outsideSpikes);
    [resultCell_{18,1},resultCell_{19,1},resultCell_{20,1}] = statPara(ISI_oSB);
else
    shortInvalid = 1;
end
%%
LB_start = burstVar.longBurst_starts;
LB_stop = burstVar.longBurst_stops;
nbLB = length(LB_start);

if length(LB_start)>2 && length(LB_stop)>2&& nbLB>2
    MLBR = nbLB/interval;
    resultCell_{21,1} = nbLB;
    resultCell_{22,1} = MLBR;
    
    ILBI_s  = spikeList(LB_start(2:end)) - spikeList(LB_stop(1:end-1));
    [resultCell_{23,1},resultCell_{24,1},resultCell_{25,1}] = statPara(ILBI_s);
    
    DLB =spikeList(LB_stop) - spikeList(LB_start);
    [resultCell_{26,1},resultCell_{27,1},resultCell_{28,1}] = statPara(DLB);
else
    longInvalid = 1;
end

if longInvalid==0 && shortInvalid==0
    [SBinside_start, SBinside_stop, nbSB_LB] = selectShinLong(SB_start, SB_stop, LB_start, LB_stop);
    [resultCell_{29,1},resultCell_{30,1},resultCell_{31,1}] = statPara(nbSB_LB);
    D_SB_LB = spikeList(SBinside_stop) - spikeList(SBinside_start);
    [resultCell_{32,1},resultCell_{33,1},resultCell_{34,1}] = statPara(D_SB_LB);
    ISBI_LB_s  = spikeList(SBinside_start(2:end)) - spikeList(SBinside_stop(1:end-1));
    [resultCell_{35,1},resultCell_{36,1},resultCell_{37,1}] = statPara(ISBI_LB_s);
    
    nbS_LB = LB_stop - LB_start;
    [resultCell_{38,1},resultCell_{39,1},resultCell_{40,1}] = statPara(nbS_LB);
    
    SBstart_oLB =burstVar.SBstart_outsideLB;
    SBstop_oLB =burstVar.SBstop_outsideLB;
    resultCell_{41,1} = length(SBstart_oLB);
    
    ISBI_oLB_s  = spikeList(SBstart_oLB(2:end)) - spikeList(SBstop_oLB(1:end-1));
    [resultCell_{42,1},resultCell_{43,1},resultCell_{44,1}] = statPara(ISBI_oLB_s);
    
    nbS_oLB =  LB_start(2:end)-1 - LB_stop(1:end-1)+1;
    resultCell_{45,1} = LB_start(1)-1 + sum(nbS_oLB) + length(spikeList)-LB_stop(end)+1;
    
    [outsideSpikes] = isolate_oSpikes(LB_start,LB_stop, spikeList);
    ISI_oLB =  diff(outsideSpikes);
    resultCell_{46,1} = resultCell_{45,1}./(interval-sum(DLB));
    [resultCell_{47,1},resultCell_{48,1},resultCell_{49,1}] = statPara(ISI_oLB);
end
%%
Cl_start = burstVar.cluster_starts;
Cl_stop = burstVar.cluster_stops;
nbCl = length(Cl_start);
if length(Cl_start)>2 && length(Cl_stop)>2 && nbCl > 2
    
    MClR = nbCl/interval;
    resultCell_{50,1} = nbCl;
    resultCell_{51,1} = MClR;
    
    IClI_s  = spikeList(Cl_start(2:end)) - spikeList(Cl_stop(1:end-1));
    [resultCell_{52,1},resultCell_{53,1},resultCell_{54,1}] = statPara(IClI_s);
    
    DCl =spikeList(Cl_stop) - spikeList(Cl_start);
    [resultCell_{55,1},resultCell_{56,1},resultCell_{57,1}] = statPara(DCl);
else
    clusterInvalid = 1;
end

if longInvalid==0 && clusterInvalid==0
    [LBinside_start, LBinside_stop, nbLB_Cl] = selectShinLong(LB_start, LB_stop, Cl_start, Cl_stop);
    [resultCell_{58,1},resultCell_{59,1},resultCell_{60,1}] = statPara(nbLB_Cl);
    D_LB_Cl = spikeList(LBinside_stop) - spikeList(LBinside_start);
    [resultCell_{61,1},resultCell_{62,1},resultCell_{63,1}] = statPara(D_LB_Cl);
    ILBI_Cl_s  = spikeList(LBinside_start(2:end)) - spikeList(LBinside_stop(1:end-1));
    [resultCell_{64,1},resultCell_{65,1},resultCell_{66,1}] = statPara(ILBI_Cl_s);
end

if longInvalid==0 && shortInvalid==0 && clusterInvalid == 0
    [SBinside_start, SBinside_stop, nbSB_Cl] = selectShinLong(SB_start, SB_stop, Cl_start, Cl_stop);
    [resultCell_{67,1},resultCell_{68,1},resultCell_{69,1}] = statPara(nbSB_Cl);
    D_SB_Cl = spikeList(SBinside_stop) - spikeList(SBinside_start);
    [resultCell_{70,1},resultCell_{71,1},resultCell_{72,1}] = statPara(D_SB_Cl);
    ISBI_Cl_s  = spikeList(SBinside_start(2:end)) - spikeList(SBinside_stop(1:end-1));
    [resultCell_{73,1},resultCell_{74,1},resultCell_{75,1}] = statPara(ISBI_Cl_s);
    
    nbS_Cl = Cl_stop - Cl_start;
    [resultCell_{76,1},resultCell_{77,1},resultCell_{78,1}] = statPara(nbS_Cl);
    
    LBstart_oCl =burstVar.LBstart_outsideCl;
    LBstop_oCl =burstVar.LBstop_outsideCl;
    resultCell_{79,1} = length(LBstart_oCl);
    ILBI_oCl_s  = spikeList(LBstart_oCl(2:end)) - spikeList(LBstop_oCl(1:end-1));
    [resultCell_{80,1},resultCell_{81,1},resultCell_{82,1}] = statPara(ILBI_oCl_s);
    
    SBstart_oCl =burstVar.SBstart_outsideCl;
    SBstop_oCl =burstVar.SBstop_outsideCl;
    resultCell_{83,1} = length(SBstart_oCl);
    ISBI_oCl_s  = spikeList(SBstart_oCl(2:end)) - spikeList(SBstop_oCl(1:end-1));
    [resultCell_{84,1},resultCell_{85,1},resultCell_{86,1}] = statPara(ISBI_oCl_s);
    
    nbS_oCl =  Cl_start(2:end)-1 - Cl_stop(1:end-1);
    resultCell_{87,1} = Cl_start(1)-1 + sum(nbS_oCl) + length(spikeList)-Cl_stop(end)+1;
    
    [outsideSpikes] = isolate_oSpikes(Cl_start,Cl_stop, spikeList);
    ISI_oCl =  diff(outsideSpikes);
    resultCell_{88,1} = resultCell_{87,1} ./(interval-sum(DCl));
    [resultCell_{89,1},resultCell_{90,1},resultCell_{91,1}] = statPara(ISI_oCl);
end
%%
Comp_start = burstVar.ClusterLB_starts;
Comp_stop = burstVar.ClusterLB_stops;
nbComp = length(Comp_start);

if ~isempty(Comp_start) && ~isempty(Comp_stop)&& ~isempty(nbComp)
    MCompR = nbComp/interval;
    resultCell_{92,1} = nbComp;
    resultCell_{93,1} = MCompR;
    
    ICompI_s  = spikeList(Comp_start(2:end)) - spikeList(Comp_stop(1:end-1));
    [resultCell_{94,1},resultCell_{95,1},resultCell_{96,1}] = statPara(ICompI_s);
    
    DComp =spikeList(Comp_stop) - spikeList(Comp_start);
    [resultCell_{97,1},resultCell_{98,1},resultCell_{99,1}] = statPara(DComp);
else
    compInvalid = 1;
end

if shortInvalid==0  && compInvalid == 0
    [SBinside_start, SBinside_stop, nbSB_Comp] = selectShinLong(SB_start, SB_stop, Comp_start, Comp_stop);
    [resultCell_{100,1},resultCell_{101,1},resultCell_{102,1}] = statPara(nbSB_Comp);
    D_SB_Comp = spikeList(SBinside_stop) - spikeList(SBinside_start);
    [resultCell_{103,1},resultCell_{104,1},resultCell_{105,1}] = statPara(D_SB_Comp);
    ISBI_Comp_s  = spikeList(SBinside_start(2:end)) - spikeList(SBinside_stop(1:end-1));
    [resultCell_{106,1},resultCell_{107,1},resultCell_{108,1}] = statPara(ISBI_Comp_s);
    
    nbS_Comp = Comp_stop - Comp_start;
    [resultCell_{109,1},resultCell_{110,1},resultCell_{111,1}] = statPara(nbS_Comp);
    
    SBstart_oComp =burstVar.SBstart_outsideCompl;
    SBstop_oComp =burstVar.SBstop_outsideCompl;
    resultCell_{112,1} = length(SBstart_oComp);
    ISBI_oComp_s  = spikeList(SBstart_oComp(2:end)) - spikeList(SBstop_oComp(1:end-1));
    [resultCell_{113,1},resultCell_{114,1},resultCell_{115,1}] = statPara(ISBI_oComp_s);
    
    nbS_oComp =  Comp_start(2:end)-1 - Comp_stop(1:end-1);
    resultCell_{116,1} = Comp_start(1)-1 + sum(nbS_oComp) + length(spikeList)-Comp_stop(end)+1;
    
    [outsideSpikes] = isolate_oSpikes(Comp_start,Comp_stop, spikeList);
    ISI_oComp =  diff(outsideSpikes);
    resultCell_{117,1} = resultCell_{116,1}./(interval-sum(DComp));
    [resultCell_{118,1},resultCell_{119,1},resultCell_{120,1}] = statPara(ISI_oComp);
end
resultCell_{121,1}= burstVar.MContrCh_SB;
resultCell_{122,1}= burstVar.CV_ContrCh_SB;
resultCell_{123,1}= burstVar.Sk_ContrCh_SB;
resultCell_{124,1}= burstVar.MContrCh_LB;
resultCell_{125,1}= burstVar.CV_ContrCh_LB;
resultCell_{126,1}= burstVar.Sk_ContrCh_LB;
resultCell_{127,1}= burstVar.MContrCh_Cl;
resultCell_{128,1}= burstVar.CV_ContrCh_Cl;
resultCell_{129,1}= burstVar.Sk_ContrCh_Cl;
resultCell_{130,1}= burstVar.MContrCh_comp;
resultCell_{131,1}= burstVar.CV_ContrCh_comp;
resultCell_{132,1}= burstVar.Sk_ContrCh_comp;
end