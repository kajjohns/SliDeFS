

GeoBounds = textread('GeoBounds.txt');

%only want SoCal segments
%ind = GeoBounds(:,3) < 37.5;
%GeoBounds = GeoBounds(ind,:);


origin = [34  -120];
%convert lat,long to x,y
x=GeoBounds(:,2:5);
llhx = fliplr(x(:,1:2));
x1 = llh2local(fliplr(llhx)', fliplr(origin))';
llhx = fliplr(x(:,3:4));
x2 = llh2local(fliplr(llhx)', fliplr(origin))';
SegEnds=[x1 x2];

SegEnds(16,1:2) = [-354.5 747];
SegEnds(1561,3:4) = [-16.21 863.8];
SegEnds(1535,1:2) = [-87.62 695.1];
SegEnds(1804,1:2) = [203.2 197.6];
SegEnds(603,1:2) = [251.4 214.8];
SegEnds(628,1:2) = [281.6 163.7];
SegEnds(633,1:2) = [291 152.2];
SegEnds(147,1:2) = [284.4 183.8];
SegEnds(317,1:2) = SegEnds(1679,3:4);
SegEnds(1673,3:4) = [318.1 61.52];
SegEnds(1488,1:2) = [-30.63 83.88];
SegEnds(1481,1:2) = [-10.5 72.4];
SegEnds(292,3:4) = [287.4 45.39];
SegEnds(1067,1:2) = [255.5 19.5];
SegEnds(378,1:2) = [303.2 11.11];
SegEnds(1012,3:4) = [289 8.412];
SegEnds(830,1:2) = SegEnds(1044,3:4);
SegEnds(337,1:2) = [367.9 22.44];
SegEnds(692,1:2) = [394.3 -62.61];
SegEnds(405,1:2) = [398.8 -84.58];
SegEnds(405,3:4) = [390.4 -96.09];
SegEnds(692,1:2) = [398 -71.1];
SegEnds(692,3:4) = [413.5 -107.2];
SegEnds(701,1:2) = [428.9 -179.8];
SegEnds(722,1:2) = [245.6 -176.4];
SegEnds(1731,1:2) = [179.5 -79.16];
SegEnds(1730,3:4) = [179.5 -79.16];
SegEnds(854,3:4) = [198.7 -12.73];
SegEnds(871,3:4) = [120.1 34.19];
SegEnds(923,3:4) = [123 46.31];
SegEnds(896,1:2) = [132.6 4.481];
SegEnds(901,1:2) = [96.02 6.09];
SegEnds(1373,1:2) = [-236.1 436.9];
SegEnds(1267,3:4) = [-223.7 404.6];
SegEnds(1149,3:4) = [-210.5 369.1];
SegEnds(1262,1:2) = [-194.7 372.4];
SegEnds(1609,1:2) = [-205.8 355.6];
SegEnds(1397,1:2) = [-171.6 349.4];
SegEnds(1145,1:2) = [-172.5 348.5];
SegEnds(1540,3:4) = [-260.5 551];
SegEnds(1539,1:2) = [-249.9 552.2];
SegEnds(1220,3:4) = [-214.2 459.2];
SegEnds(1467,1:2) = [-206.1 525.9];
SegEnds(1193,1:2) = SegEnds(1158,1:2);
SegEnds(1188,3:4) = [-185.7 435.3];
SegEnds(1234,1:2) = [-149.4 410];
SegEnds(1237,3:4) = [-181.6 454.4];
SegEnds(1172,1:2) = [-177.5 402.9];
SegEnds(1207,1:2) = [-148.6 360.8];
SegEnds(1163,1:2) = SegEnds(1171,3:4);
SegEnds(1402,3:4) = [-129 321.1];
SegEnds(1619,3:4) = [-133.5 312];
SegEnds(588,3:4) = SegEnds(613,1:2);
SegEnds(1,1:2) = [-101.4 340.8];
SegEnds(613,3:4) = SegEnds(589,1:2);
SegEnds(589,3:4) = SegEnds(590,1:2);
SegEnds(591,3:4) = SegEnds(595,1:2);
SegEnds(1605,3:4) = [-107.3 186.1];
SegEnds(223,3:4) = [-85.1 148.2];
SegEnds(1144,3:4) = [-83.4 137.8];
SegEnds(223,1:2) = [-87.99 149.9];
SegEnds(1300,3:4) = SegEnds(1456,3:4);
SegEnds(230,3:4) = [-42.69 123.2];
SegEnds(1594,1:2) = SegEnds(1515,1:2);
SegEnds(1449,3:4) = SegEnds(1488,1:2);
SegEnds(1481,1:2) = [-13.04 73.31];
SegEnds(1781,3:4) = [7.343 61.7];
SegEnds(1774,3:4) = [10.78 62.73];
SegEnds(1786,3:4) = [59.2 76.23];
SegEnds(1055,3:4) = [99.43 90.56];
SegEnds(1049,1:2) = [99.43 90.56];
SegEnds(815,1:2) = [88.76 89.61];
SegEnds(923,3:4) = [121.7 47.28];
SegEnds(849,3:4) = [81.12 48.02];
SegEnds(921,3:4) = [131 48.34];
SegEnds(506,1:2) = SegEnds(871,3:4);
SegEnds(612,3:4) = [158.5 12.56];
SegEnds(620,3:4) = [175.3 15.66];
SegEnds(621,1:2) = [175.3 15.66];
SegEnds(479,1:2) = [171.1 15.61];
SegEnds(954,3:4) = [168.2 -.6957];
SegEnds(924,3:4) = [145.8 -20.23];
SegEnds(460,3:4) = [145.8 -15.88];
SegEnds(461,1:2) = [145.8 -15.88];
SegEnds(584,1:2) = [-192 557.2];
SegEnds(1839,3:4) = [253.7 251.8];
SegEnds(1008,3:4) = [312.1366 7];
SegEnds(1009,1:2) = [312.1366 7];
SegEnds(1008,1:2) = [316.8 5.5];
SegEnds(558,3:4) = SegEnds(497,3:4);
SegEnds(614,1:2) = SegEnds(496,3:4);
SegEnds(877,2) = SegEnds(877,2)-1;
SegEnds(877,4) = SegEnds(877,4)-1;
SegEnds(878,2) = SegEnds(878,2)-1;
SegEnds(878,4) = SegEnds(878,4)-1;
SegEnds(879,2) = SegEnds(879,2)-1;



%only want SoCal segments
ind = GeoBounds(:,3) < 37.5;
GeoBounds = GeoBounds(ind,:);
SegEnds = SegEnds(ind,:);

dip = GeoBounds(:,9);
dip(dip==0)=90;


%calculate lengths of segments
SegLength=sqrt((SegEnds(:,1)-SegEnds(:,3)).^2+(SegEnds(:,2)-SegEnds(:,4)).^2);
%calculate strike of segments
angle=atan2(SegEnds(:,4)-SegEnds(:,2),SegEnds(:,3)-SegEnds(:,1));
strike=90-angle*180/pi;
%calculate center of segments
centers=[(SegEnds(:,1)+SegEnds(:,3))./2 (SegEnds(:,2)+SegEnds(:,4))./2];
SegEndsllh=x;
centersllh=[(SegEndsllh(:,1)+SegEndsllh(:,3))./2 (SegEndsllh(:,2)+SegEndsllh(:,4))./2];
%get vector normal to segments
SegNorm=[-sin(angle) cos(angle)];
%get vector parallel to segments
SegStrike=[cos(angle) sin(angle)];


seg_strike = strike;

%locking depth
H=toprow_depth;
nve=1;


%faults=[length, width, *depth, dip, strike(degrees), *north offset, *east offset]
% *depth to top edge, north and east offsets refer to location of center of top edge
widths=H./sin(dip*pi/180);
temp=abs(widths).*cos(dip*pi/180);
xoffset=-temp.*cos(pi/2+angle);
yoffset=-temp.*sin(pi/2+angle);
faults=[SegLength widths H*ones(size(SegLength)) dip strike .5*(SegEnds(:,1)+SegEnds(:,3))+xoffset .5*(SegEnds(:,2)+SegEnds(:,4))+yoffset];


nhe = ceil(SegLength/pL);  %number horizontal elements 


ub_Segs = GeoBounds(:,7);
lb_Segs = GeoBounds(:,6);
rake = GeoBounds(:,8);


%something wrong with seg 1034
ub_Segs(1034)=18;
lb_Segs(1034)=10;


pm_top = [];

ub_top_patches=[];
lb_top_patches=[];
rake_top_patches=[];
dip_top_patches=[];
strike_top_patches=[];
centers_top_patches=[];

pm_top_seg_num = [];

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
rake_top_patches = [rake_top_patches; rake(k)*ones(nhe(k)*nve,1)];
dip_top_patches = [dip_top_patches; dip(k)*ones(nhe(k)*nve,1)];
strike_top_patches = [strike_top_patches; strike(k)*ones(nhe(k)*nve,1)];
centers_top_patches = [centers_top_patches; [pf(:,6)-xoffset(k)/nve pf(:,7)-yoffset(k)/nve]];

pm_top_seg_num = [pm_top_seg_num; k*ones(size(pf,1),1)];

end




PatchEnds_top = [centers_top_patches(:,1)+.5*pm_top(:,1).*cos( (90-strike_top_patches)*pi/180) centers_top_patches(:,2)+.5*pm_top(:,1).*sin( (90-strike_top_patches)*pi/180) centers_top_patches(:,1)-.5*pm_top(:,1).*cos( (90-strike_top_patches)*pi/180) centers_top_patches(:,2)-.5*pm_top(:,1).*sin( (90-strike_top_patches)*pi/180)];





 