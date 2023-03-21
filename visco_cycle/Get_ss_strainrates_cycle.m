function [Gss, Gds] = Get_ss_strainrates_cycle(Segs,dips,Ld,Cd,xystats,...
                      H1,get_SS,get_DS)

                 

if nargin<14
    Nterms=[];
end



xloc = [xystats';zeros(1,size(xystats,1))];
get_SS = logical(get_SS);
get_DS = logical(get_DS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%segment geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 centers = [(Segs(:,1)+Segs(:,2))/2 (Segs(:,3)+Segs(:,4))/2];
%calculate strike of segments
angle=atan2(Segs(:,4)-Segs(:,3),Segs(:,2)-Segs(:,1));
strike=90-angle*180/pi;
SegLength = sqrt( (Segs(:,4)-Segs(:,3)).^2 + (Segs(:,2)-Segs(:,1)).^2 );

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



for k=1:size(centers,1)  %loop over number of segments
    
        W = widths(k);  %width of fault
        
        
        %strike-slip contribution 
        m_ss=[SegLength(k) W H1 dips(k) strike(k) bottcenters(k,1) bottcenters(k,2) 1 0 0];
        %dip-slip contribution 
        m_ds=[SegLength(k) W H1 dips(k) strike(k) bottcenters(k,1) bottcenters(k,2) 0 1 0];
     
      
           
        
            if get_SS(k)  %only compute if requested -- it is expensive to compute and not use 
            
                
                
                %compute velocities at four points around the observation
                %points and use difference to compute strain rates
                dx = 1; dy = 1;  %1 km
                
                xloc_plusx = xloc; xloc_plusx(1,:) = xloc_plusx(1,:)+dx;
                xloc_minusx = xloc; xloc_minusx(1,:) = xloc_minusx(1,:)-dx;
                
                xloc_plusy = xloc; xloc_plusy(2,:) = xloc_plusy(2,:)+dy;
                xloc_minusy = xloc; xloc_minusy(2,:) = xloc_minusy(2,:)-dy;
              
                Uss_plusx = LayeredBasis_viscous(m_ss,xloc_plusx,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uss_minusx = LayeredBasis_viscous(m_ss,xloc_minusx,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uss_plusy = LayeredBasis_viscous(m_ss,xloc_plusy,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uss_minusy = LayeredBasis_viscous(m_ss,xloc_minusy,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
            
                
                Dx =  (Uss_plusx - Uss_minusx)/(2*dx);
                Dy =  (Uss_plusy - Uss_minusy)/(2*dy);
                
                Exx_ss = Dx(1,:)';
                Exy_ss = .5*(Dx(2,:)+Dy(1,:))';
                Eyy_ss = Dy(2,:)';
                
                [Ue,Dss,S,flag]=disloc3d(m_ss',xloc,1,.25); %elastic part
                
                Exx_ss = Exx_ss + Dss(1,:)';
                Exy_ss = Exy_ss + .5*(Dss(2,:)+Dss(4,:))';
                Eyy_ss = Eyy_ss + Dss(5,:)';
                
            else
                
                Exx_ss = zeros(size(xloc,2),1);
                Exy_ss = zeros(size(xloc,2),1);
                Eyy_ss = zeros(size(xloc,2),1);
                
            end
        
            if get_DS(k)  %only compute if requested -- it is expensive to compute and not use 
                
                %compute velocities at four points around the observation
                %points and use difference to compute strain rates
                dx = 1; dy = 1;  %1 km
                
                xloc_plusx = xloc; xloc_plusx(1,:) = xloc_plusx(1,:)+dx;
                xloc_minusx = xloc; xloc_minusx(1,:) = xloc_minusx(1,:)-dx;
                
                xloc_plusy = xloc; xloc_plusy(2,:) = xloc_plusy(2,:)+dy;
                xloc_minusy = xloc; xloc_minusy(2,:) = xloc_minusy(2,:)-dy;
              
                Uds_plusx = LayeredBasis_viscous(m_ds,xloc_plusx,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uds_minusx = LayeredBasis_viscous(m_ds,xloc_minusx,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uds_plusy = LayeredBasis_viscous(m_ds,xloc_plusy,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
                Uds_minusy = LayeredBasis_viscous(m_ds,xloc_minusy,H1,1,1,10^6,1);  %note positive strike-slip is left-lateral
            
                
                Dx =  (Uds_plusx - Uds_minusx)/(2*dx);
                Dy =  (Uds_plusy - Uds_minusy)/(2*dy);
                
                Exx_ds = Dx(1,:)';
                Exy_ds = .5*(Dx(2,:)+Dy(1,:))';
                Eyy_ds = Dy(2,:)';
                
                [Ue,Dds,S,flag]=disloc3d(m_ds',xloc,1,.25); %elastic part
                
                Exx_ds = Exx_ds + Dds(1,:)';
                Exy_ds = Exy_ds + .5*(Dds(2,:)+Dds(4,:))';
                Eyy_ds = Eyy_ds + Dds(5,:)';
            
            else
                
                Exx_ds = zeros(size(xloc,2),1);
                Exy_ds = zeros(size(xloc,2),1);
                Eyy_ds = zeros(size(xloc,2),1);
                
            end
        
        
     
  
     
    Gss.Exx(:,k) = Exx_ss;
    Gss.Exy(:,k) = Exy_ss;
    Gss.Eyy(:,k) = Eyy_ss;
  
    Gds.Exx(:,k) = Exx_ds;
    Gds.Exy(:,k) = Exy_ds;
    Gds.Eyy(:,k) = Eyy_ds;
  
     
end %for k=1:size(centers,1)



Gss.Exx = sum(Gss.Exx,2);
Gss.Exy = sum(Gss.Exy,2);
Gss.Eyy = sum(Gss.Eyy,2);

Gds.Exx = sum(Gds.Exx,2);
Gds.Exy = sum(Gds.Exy,2);
Gds.Eyy = sum(Gds.Eyy,2);



