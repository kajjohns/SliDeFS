
%% This script builds "edges" of fault patch model along fault segments for use in constructin the triangular mesh
% Note, the actual fault patch model is constructed with a diffferent
% script


%calculate lengths of segments
SegLength=sqrt((SegEnds(:,1)-SegEnds(:,3)).^2+(SegEnds(:,2)-SegEnds(:,4)).^2);
%calculate strike of segments
angle=atan2(SegEnds(:,4)-SegEnds(:,2),SegEnds(:,3)-SegEnds(:,1));
strike=90-angle*180/pi;
%calculate center of segments
centers=[(SegEnds(:,1)+SegEnds(:,3))./2 (SegEnds(:,2)+SegEnds(:,4))./2];
%get vector normal to segments
SegNorm=[-sin(angle) cos(angle)];
%get vector parallel to segments
SegStrike=[cos(angle) sin(angle)];


seg_strike = strike;


%faults=[length, width, *depth, dip, strike(degrees), *north offset, *east offset]
% *depth to top edge, north and east offsets refer to location of center of top edge
nve=1;
Ones = 90*ones(size(SegEnds,1),1);
faults=[SegLength Ones Ones 90*Ones strike .5*(SegEnds(:,1)+SegEnds(:,3)) .5*(SegEnds(:,2)+SegEnds(:,4))];


	
nhe = ceil(SegLength/pL);  %number horizontal elements 


edge = [];
node =[];
for k=1:size(faults,1)

%specify components of slip to be calculate ([strike-slip,dip-slip,opening]) -- e.g. [0 1 0] means dip slip only
dis_geom  = [faults(k,:), [1 1 0]];
%dis_geom = movefault(dis_geom1);  % move the fault so that the coordinates of the midpoint refer to the
											 % fault bottom as in Okada

	
%% Create slip patches
pm = patchfault(dis_geom(1,1:7),nhe(k),nve);


segends1 = [pm(:,6)+ pm(:,1)/2.*cos((90-pm(:,5))*pi/180) pm(:,7)+ pm(:,1)/2.*sin((90-pm(:,5))*pi/180)];
segends2 = [pm(:,6)- pm(:,1)/2.*cos((90-pm(:,5))*pi/180) pm(:,7)- pm(:,1)/2.*sin((90-pm(:,5))*pi/180)];

%segends = [segends1;segends2];
segends = zeros(size(segends1,1)+size(segends2,1),2);
segends(1:2:end) = segends1;
segends(2:2:end) = segends2;

node = [node; segends];

edge1 = [(1:2:size(segends,1))' (2:2:size(segends,1))'];

if k>1
    edge = [edge; edge1+edge(end,end)];
else
    edge = edge1;
end



end






 