
%unwrap all fault info and arrange in matrices/vectors
SegEnds = [];
Seg_FaultID = [];
Seg_miniID = [];
Seg_Dips= [];
Seg_DipDir= [];
Seg_Rakes = [];
Seg_LowDepth = [];
Seg_UpDepth = [];
Seg_GeoPrefRate = [];
Seg_GeoLowRate = [];
Seg_GeoHighRate = [];
Seg_TimeSinceEq = [];
Seg_RecurrenceTime = [];

std_slip_rate = []; 

for k=1:length(FaultInfo)

    N = size(FaultInfo{k}.SegEnds,1);

    SegEnds = [SegEnds;FaultInfo{k}.SegEnds];
    Seg_FaultID = [Seg_FaultID; repmat(FaultInfo{k}.FaultID,N,1)];
    

    Seg_Dips = [Seg_Dips; repmat(FaultInfo{k}.Dip,N,1)];
    Seg_Rakes = [Seg_Rakes; repmat(FaultInfo{k}.Rake,N,1)];
    
    Seg_GeoPrefRate = [Seg_GeoPrefRate; repmat(FaultInfo{k}.PrefSlipRate,N,1)];
    Seg_GeoLowRate = [Seg_GeoLowRate; repmat(FaultInfo{k}.LBSlipRate,N,1)];
    Seg_GeoHighRate = [Seg_GeoHighRate; repmat(FaultInfo{k}.UBSlipRate,N,1)];

    std_slip_rate = [std_slip_rate; repmat(FaultInfo{k}.STDSlipRate,N,1)];
    
    Seg_RecurrenceTime = [Seg_RecurrenceTime; repmat(FaultInfo{k}.RecurrenceTime,N,1)];
    Seg_TimeSinceEq = [Seg_TimeSinceEq; repmat(FaultInfo{k}.TimeSinceEq,N,1)];

   
end


pref_slip_rate = Seg_GeoPrefRate;
ub_slip_rate = Seg_GeoHighRate;
lb_slip_rate = Seg_GeoLowRate;

%set empty matrices to nans
if isempty(pref_slip_rate)
    pref_slip_rate = nan(size(Seg_Dips));
end
if isempty(ub_slip_rate)
    ub_slip_rate = nan(size(Seg_Dips));
end
if isempty(lb_slip_rate)
    lb_slip_rate = nan(size(Seg_Dips));
end
if isempty(std_slip_rate)
    std_slip_rate = nan(size(Seg_Dips));
end


