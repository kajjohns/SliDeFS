nve=nve_total;
H = DLock_max;


%faults=[length, width, *depth, dip, strike(degrees), *north offset, *east offset]
% *depth to top edge, north and east offsets refer to location of center of top edge
widths=H./sin(dip*pi/180);
temp=abs(widths).*cos(dip*pi/180);
xoffset=-temp.*cos(pi/2+angle);
yoffset=-temp.*sin(pi/2+angle);
faults=[SegLength widths H dip strike .5*(SegEnds(:,1)+SegEnds(:,3))+xoffset .5*(SegEnds(:,2)+SegEnds(:,4))+yoffset];


nhe = ceil(SegLength/pL);  %number horizontal elements 

pm = [];

ub_patches=[];
lb_patches=[];
pref_patches=[];
std_patches=[];
rake_patches=[];
dip_patches=[];
strike_patches=[];
centers_patches=[];

pm_seg_num = [];
pm_seg_id = [];

for k=1:size(faults,1)

%specify components of slip to be calculate ([strike-slip,dip-slip,opening]) -- e.g. [0 1 0] means dip slip only
dis_geom  = [faults(k,:), [1 1 0]];
%dis_geom = movefault(dis_geom1);  % move the fault so that the coordinates of the midpoint refer to the
											 % fault bottom as in Okada

	
%% Create slip patches
pf = patchfault(dis_geom(1,1:7),nhe(k),nve);

%toss out top row
pf(1:nve:end,:) = [];


pm = [pm; pf];

ub_patches = [ub_patches; ub_Segs(k)*ones(nhe(k)*(nve-1),1)];
lb_patches = [lb_patches; lb_Segs(k)*ones(nhe(k)*(nve-1),1)];
pref_patches = [pref_patches; pref_Segs(k)*ones(nhe(k)*(nve-1),1)];
std_patches = [std_patches; std_Segs(k)*ones(nhe(k)*(nve-1),1)];
rake_patches = [rake_patches; rake(k)*ones(nhe(k)*(nve-1),1)];
dip_patches = [dip_patches; dip(k)*ones(nhe(k)*(nve-1),1)];
strike_patches = [strike_patches; strike(k)*ones(nhe(k)*(nve-1),1)];
centers_patches = [centers_patches; [pf(:,6)-xoffset(k)/(nve-1) pf(:,7)-yoffset(k)/(nve-1)]];

pm_seg_num = [pm_seg_num; k*ones(size(pf,1),1)];
pm_seg_id = [pm_seg_id; segid(k)*ones(size(pf,1),1)];

end




PatchEnds = [centers_patches(:,1)+.5*pm(:,1).*cos( (90-strike_patches)*pi/180) centers_patches(:,2)+.5*pm(:,1).*sin( (90-strike_patches)*pi/180) centers_patches(:,1)-.5*pm(:,1).*cos( (90-strike_patches)*pi/180) centers_patches(:,2)-.5*pm(:,1).*sin( (90-strike_patches)*pi/180)];





 