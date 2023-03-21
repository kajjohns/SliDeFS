function [Gss, Gds] = Get_steady_vels(Segs,dips,xystats,xy_visco,H1,get_SS,get_DS)

                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  [Gss, Gds] = Get_interseismic_vels_cycle(Segs,dips,Ld,Cd,xystats,...
%                      H1,H2,tR1,tR2,T,timeEq,get_SS,get_DS)build_layered_cycle_greens.m
%
%  This function builds the greens functions response for periodic
%  earthquakes on rectangular faults in an elastic 
%  layer overlying a (Maxwell) viscoelastic layer over a (Maxwell) viscoelastic 
%  half space. Outputs are velocites at specified time and specified
%  locations. This computes the response to steady backslip on the fault
%  (to cancel the long-term rate and lock the fault) plus periodic
%  earthquakes.
%
%  INPUTS:  
%
%  H1 =  thickness of elastic layer  (km)
%  H2 =  depth to bottom of viscoelastic layer -- km -- (same as top of half-space)
%  tR1 = relaxation time (2*eta/mu, eta=viscosity, mu=shear modulus) of viscoelastic layer -- years
%  tR2 = relaxation time of viscoelastic halfspace -- years
%  T = vector of earthquake recurrence times (years) -- one for each
%      rectangular section  (NOTE:  NaN indicates no earthquakes and no viscoelastic effect --
%      compute only elastic response to back slip)
%  timeEq  =  vector of times since last earthquake (years) one for each
%      rectangular segment NOTE:  NaN indicates no earthquakes and no viscoelastic effect --
%      compute only elastic response to back slip)
%  Ends =  matrix of Endpoints of surface trace of rectangular segments (km)
%          each row is a segment:  x1  x2  y1  y2  (where x is East, y is
%          north and 1 and 2 refer to first and second endpoint)
%  dips =  vector of segment dips (degrees) using right-hand rule (positive
%          dip down to right when facing from endpoint1 to endpoint2)
%  Ld = vector of locking depths (creep at constant stress below)
%  Cd = vector of creepding depth (crep at steady rate above to surface)
%          Cd=0 if no surface creep; of course Ld>=Cd
%  xystats = Nx2 matrix of (x,y) coordinates of N observation points
%  get_SS = vector of ones and zeros, specifying whether or not to
%           calculate contribution from strike-component of slip
%           (1=yes, 0=no) -- if 'no', returns NaN
%  get_DS = vector of ones and zeros, specifying whether or not to
%           calculate contribution from dip-component of slip
%           (1=yes, 0=no) -- if 'no', returns NaN
%  Nterms = number of terms in numerical Hankel transform (default=100 if not
%           specified)
%
%  OUTPUTS:
%
%  Gss = 3NxM vector of velocities at N observation points for M
%        rectangular sources (unit strike-slip rate imposed on each source)
%        Gss(1:end/3) = east component
%        Gss(1+end/3:2*end/3) = north component
%        Gss(1+2*end/3:end) = vertical component
%
%  Gds = same as above for unit dip-slip rate imposed on each source
%
%  Kaj Johnson, Indiana University, 2007-2012, last modified Jan. 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin<14
    Nterms=[];
end



xloc = [xystats';zeros(1,size(xystats,1))];
xloc_visco_grid = [xy_visco';zeros(1,size(xy_visco,1))];

get_SS = logical(get_SS);
get_DS = logical(get_DS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%segment geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 centers = [(Segs(:,1)+Segs(:,3))/2 (Segs(:,2)+Segs(:,4))/2];
%calculate strike of segments
angle=atan2(Segs(:,4)-Segs(:,2),Segs(:,3)-Segs(:,1));
strike=90-angle*180/pi;
SegLength = sqrt( (Segs(:,4)-Segs(:,2)).^2 + (Segs(:,3)-Segs(:,1)).^2 );

%get center of bottom segment

%this next part is only necessary for dipping faults

%get center of bottom segment
widths=H1./sin(dips*pi/180);
temp=abs(widths).*cos(dips*pi/180);
xoffset=-temp.*cos(pi/2+angle);
yoffset=-temp.*sin(pi/2+angle);
bottcenters=centers+[xoffset yoffset];



%get vector normal to plane
StrikeVec=[cos(angle) sin(angle) 0*angle];
DipVec=[bottcenters-centers -H1*ones(size(bottcenters,1),1)];
NormalVec=cross(StrikeVec,DipVec,2);
for k=1:size(NormalVec,1)
    NormalVec(k,:)=NormalVec(k,:)/norm(NormalVec(k,:));
    DipVec(k,:)=DipVec(k,:)/norm(DipVec(k,:));
end



Gss.east=zeros(size(xystats,1),size(centers,1));
Gss.north=zeros(size(xystats,1),size(centers,1));
Gss.up=zeros(size(xystats,1),size(centers,1));

Gds.east=zeros(size(xystats,1),size(centers,1));
Gds.north=zeros(size(xystats,1),size(centers,1));
Gds.up=zeros(size(xystats,1),size(centers,1));

for k=1:size(centers,1)  %loop over number of segments
    
        W = widths(k);  %width of fault
        
        
        %strike-slip contribution 
        m_ss=[SegLength(k) W H1 dips(k) strike(k) bottcenters(k,1) bottcenters(k,2) 1 0 0];
        %dip-slip contribution 
        m_ds=[SegLength(k) W H1 dips(k) strike(k) bottcenters(k,1) bottcenters(k,2) 0 1 0];
     
      
           
        
            if get_SS(k)  %only compute if requested -- it is expensive to compute and not use 
            
                [Ue,D,S,flag]=disloc3d(m_ss',xloc_visco_grid,1,.25); %elastic part
                U = Ue + LayeredBasis_viscous(m_ss,xloc_visco_grid,H1,1,1,10^6,1);  
            
                %interpolate to desired grid
                Uss(1,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(1,:),xloc(1,:),xloc(2,:));
                Uss(2,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(2,:),xloc(1,:),xloc(2,:));
                Uss(3,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(3,:),xloc(1,:),xloc(2,:));
                Uss(isnan(Uss))=0;


            else
                Uss = zeros(3,size(xloc,2));
            end
        
            if get_DS(k)  %only compute if requested -- it is expensive to compute and not use 
                [Ue,D,S,flag]=disloc3d(m_ds',xloc_visco_grid,1,.25); %elastic part
                U = Ue + LayeredBasis_viscous(m_ds,xloc_visco_grid,H1,1,1,10^6,1);

                 %interpolate to desired grid
                Uds(1,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(1,:),xloc(1,:),xloc(2,:));
                Uds(2,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(2,:),xloc(1,:),xloc(2,:));
                Uds(3,:) = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),U(3,:),xloc(1,:),xloc(2,:));
                Uds(isnan(Uss))=0;

            else
                Uds = zeros(3,size(xloc,2));
            end
        
        
     

     
    Gss.east(:,k) = Uss(1,:)';
    Gss.north(:,k) = Uss(2,:)';    
    Gss.up(:,k) = Uss(3,:)';
    
    Gds.east(:,k) = Uds(1,:)';
    Gds.north(:,k) = Uds(2,:)';    
    Gds.up(:,k) = Uds(3,:)';

   
   
end %for k=1:size(centers,1)


%Gss = sum(Gss,2);
%Gds = sum(Gds,2);




