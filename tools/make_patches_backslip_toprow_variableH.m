

dip = Dips;
dip(dip==0)=90;



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



%locking depth
H = DLock_max/nve_total;
nve=1;


%faults=[length, width, *depth, dip, strike(degrees), *north offset, *east offset]
% *depth to top edge, north and east offsets refer to location of center of top edge
widths=H./sin(dip*pi/180);
temp=abs(widths).*cos(dip*pi/180);
xoffset=-temp.*cos(pi/2+angle);
yoffset=-temp.*sin(pi/2+angle);
faults=[SegLength widths H dip strike .5*(SegEnds(:,1)+SegEnds(:,3))+xoffset .5*(SegEnds(:,2)+SegEnds(:,4))+yoffset];


nhe = ceil(SegLength/pL);  %number horizontal elements 




ub_Segs = ub_slip_rate;
lb_Segs =lb_slip_rate;
pref_Segs = pref_slip_rate;
std_Segs = std_slip_rate;
rake = Seg_Rakes;
segid = Seg_FaultID;


pm_top = [];

ub_top_patches=[];
lb_top_patches=[];
pref_top_patches=[];
std_top_patches=[];
rake_top_patches=[];
dip_top_patches=[];
strike_top_patches=[];
centers_top_patches=[];
pm_top_seg_num = [];
pm_top_seg_id = [];

for k=1:size(faults,1)

%specify components of slip to be calculate ([strike-slip,dip-slip,opening]) -- e.g. [0 1 0] means dip slip only
dis_geom  = [faults(k,:), [1 1 0]];
%dis_geom = movefault(dis_geom1);  % move the fault so that the coordinates of the midpoint refer to the
											 % fault bottom as in Okada

	
%% Create slip patches
pf = patchfault(dis_geom(1,1:7),nhe(k),nve);
pm_top = [pm_top; pf];

ub_top_patches = [ub_top_patches; ub_Segs(k)*ones(nhe(k)*nve,1)];
lb_top_patches = [lb_top_patches; lb_Segs(k)*ones(nhe(k)*nve,1)];
pref_top_patches = [pref_top_patches; pref_Segs(k)*ones(nhe(k)*nve,1)];
std_top_patches = [std_top_patches; std_Segs(k)*ones(nhe(k)*nve,1)];
rake_top_patches = [rake_top_patches; rake(k)*ones(nhe(k)*nve,1)];
dip_top_patches = [dip_top_patches; dip(k)*ones(nhe(k)*nve,1)];
strike_top_patches = [strike_top_patches; strike(k)*ones(nhe(k)*nve,1)];
centers_top_patches = [centers_top_patches; [pf(:,6)-xoffset(k)/nve pf(:,7)-yoffset(k)/nve]];

pm_top_seg_num = [pm_top_seg_num; k*ones(size(pf,1),1)];
pm_top_seg_id = [pm_top_seg_id; segid(k)*ones(size(pf,1),1)];

end




PatchEnds_top = [centers_top_patches(:,1)+.5*pm_top(:,1).*cos( (90-strike_top_patches)*pi/180) centers_top_patches(:,2)+.5*pm_top(:,1).*sin( (90-strike_top_patches)*pi/180) centers_top_patches(:,1)-.5*pm_top(:,1).*cos( (90-strike_top_patches)*pi/180) centers_top_patches(:,2)-.5*pm_top(:,1).*sin( (90-strike_top_patches)*pi/180)];





 