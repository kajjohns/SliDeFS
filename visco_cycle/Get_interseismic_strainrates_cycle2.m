function [Gss, Gds,Gss_e, Gds_e] = Get_interseismic_strainrates_cycle(Segs,dips,Ld,Cd,xystats,xs,ys,...
                      H1,H2,tR1,tR2,T,timeEq,get_SS,get_DS,Nterms)

                 

%compute viscoelastic response on a radially decreasgin grid size
r = logspace(0,3,15);
dth = 2*pi/15; theta = 0:dth:2*pi-dth;
[R,TH] = meshgrid(r,theta);
xgrid = R.*cos(TH);  xgrid=xgrid(:);
ygrid = R.*sin(TH);  ygrid=ygrid(:);
                  
if nargin<14
    Nterms=[];
end



xloc = [xystats';zeros(1,size(xystats,1))];
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

%get center of top segment
% widths=Cd./sin(dips*pi/180);
% temp=abs(widths).*cos(dips*pi/180);
% xoffset=-temp.*cos(pi/2+angle);
% yoffset=-temp.*sin(pi/2+angle);
% topcenters=centers+[xoffset yoffset];


%get center of bottom of locked section
widths=Ld./sin(dips*pi/180);
temp=abs(widths).*cos(dips*pi/180);
xoffset=-temp.*cos(pi/2+angle);
yoffset=-temp.*sin(pi/2+angle);
bottcenters_locked=centers+[xoffset yoffset];



%get vector normal to plane
StrikeVec=[cos(angle) sin(angle) 0*angle];
DipVec=[bottcenters-centers -H1*ones(size(bottcenters,1),1)];
NormalVec=cross(StrikeVec,DipVec,2);
for k=1:size(NormalVec,1)
    NormalVec(k,:)=NormalVec(k,:)/norm(NormalVec(k,:));
    DipVec(k,:)=DipVec(k,:)/norm(DipVec(k,:));
end



%Gss=zeros(3*size(xystats,1),size(centers,1));
%Gds=zeros(3*size(xystats,1),size(centers,1));


for k=1:size(centers,1)  %loop over number of segments
    
        W = (Ld(k)-Cd(k))./sin(dips(k)*pi/180);  %width of locked part
        
        %strike-slip contribution 
        m_ss=[SegLength(k) W Ld(k)+0.001 dips(k) strike(k) bottcenters_locked(k,1) bottcenters_locked(k,2) 1 0 0];
        %dip-slip contribution 
        m_ds=[SegLength(k) W Ld(k)+0.001 dips(k) strike(k) bottcenters_locked(k,1) bottcenters_locked(k,2) 0 1 0];
     
      

    if isnan(T(k))  %elastic part only
           
             if get_SS(k)  %only compute if requested -- it is expensive to compute and not use 
                [Uss,Dss] = backslip_elastic3D(m_ss,xloc);
                
                Exx_ss = Dss(1,:)';
                Exy_ss = .5*(Dss(2,:)+Dss(4,:))';
                Eyy_ss = Dss(5,:)';
                
                Exx_ss_elastic = Exx_ss;
                Exy_ss_elastic = Exy_ss;
                Eyy_ss_elastic = Eyy_ss;
                
                
             else
                 
                Exx_ss = zeros(size(xloc,2),1);
                Exy_ss = zeros(size(xloc,2),1);
                Eyy_ss = zeros(size(xloc,2),1);
                
                Exx_ss_elastic = Exx_ss;
                Exy_ss_elastic = Exy_ss;
                Eyy_ss_elastic = Eyy_ss;
                
                
             end

             if get_DS(k)  %only compute if requested -- it is expensive to compute and not use 
                  [Uds,Dds] = backslip_elastic3D(m_ds,xloc);
                  
                  Exx_ds = Dds(1,:)';
                  Exy_ds = .5*(Dds(2,:)+Dds(4,:))';
                  Eyy_ds = Dds(5,:)';
                  
                Exx_ds_elastic = Exx_ds;
                Exy_ds_elastic = Exy_ds;
                Eyy_ds_elastic = Eyy_ds;
                
             else
                 
                  Exx_ds = zeros(size(xloc,2),1);
                  Exy_ds = zeros(size(xloc,2),1);
                  Eyy_ds = zeros(size(xloc,2),1);
                  
                  Exx_ds_elastic = Exx_ds;
                  Exy_ds_elastic = Exy_ds;
                  Eyy_ds_elastic = Eyy_ds;
                
                  
             end
             
    else  %compute viscoelastic and elastic parts
        
            if get_SS(k)  %only compute if requested -- it is expensive to compute and not use 
                %note calculation below has no elastic backslip. This needs
                %to be added back on
                
                %compute velocities at four points around the observation
                %points and use difference to compute strain rates
                dx = 1; dy = 1;  %1 km
                
                 keyboard
                
                Uss_plusx=backslip_cycle3D_layered_noelastic(m_ss,xloc_plusx,H1,H2,1,1,timeEq(k),tR1,tR2,T(k),Nterms);
                Uss_minusx=backslip_cycle3D_layered_noelastic(m_ss,xloc_minusx,H1,H2,1,1,timeEq(k),tR1,tR2,T(k),Nterms);
                Uss_plusy=backslip_cycle3D_layered_noelastic(m_ss,xloc_plusy,H1,H2,1,1,timeEq(k),tR1,tR2,T(k),Nterms);
                Uss_minusy=backslip_cycle3D_layered_noelastic(m_ss,xloc_minusy,H1,H2,1,1,timeEq(k),tR1,tR2,T(k),Nterms);
                
         
                Dx =  (Uss_plusx - Uss_minusx)/(2*dx);
                Dy =  (Uss_plusy - Uss_minusy)/(2*dy);
                
                Exx_ss = Dx(1,:)';
                Exy_ss = .5*(Dx(2,:)+Dy(1,:))';
                Eyy_ss = Dy(2,:)';
                
                %add elastic backslip contribution
                 [Uss,Dss] = backslip_elastic3D(m_ss,xloc);
             
                   Exx_ss = Exx_ss + Dss(1,:)';
                   Exy_ss = Exy_ss + .5*(Dss(2,:)+Dss(4,:))';
                   Eyy_ss = Eyy_ss + Dss(5,:)';
             
                   Exx_ss_elastic = Dss(1,:)';
                   Exy_ss_elastic = .5*(Dss(2,:)+Dss(4,:))';
                   Eyy_ss_elastic = Dss(5,:)';
            else
                
                Exx_ss = zeros(size(xloc,2),1);
                Exy_ss = zeros(size(xloc,2),1);
                Eyy_ss = zeros(size(xloc,2),1);
                
                Exx_ss_elastic = Exx_ss;
                Exy_ss_elastic = Exy_ss;
                Eyy_ss_elastic = Eyy_ss;
                
            end
        
            if get_DS(k)  %only compute if requested -- it is expensive to compute and not use 
                
                 
                %compute velocities at four points around the observation
                %points and use difference to compute strain rates
                dx = 1; dy = 1;  %1 km
                
                
                keyboard
                
                xloc_radial_grid = [xgrid'+centers(k,1);ygrid'+centers(k,2);zeros(1,length(xgrid))];
                
                Uds_radial_grid=backslip_cycle3D_layered_noelastic(m_ds,xloc_radial_grid,H1,H2,1,1,timeEq(k),tR1,tR2,T(k),Nterms);
                
                Uds(1,:) = griddata(xloc_radial_grid(1,:),xloc_radial_grid(2,:),Uds_radial_grid(1,:),xloc(1,:),xloc(2,:));
                Uds(2,:) = griddata(xloc_radial_grid(1,:),xloc_radial_grid(2,:),Uds_radial_grid(2,:),xloc(1,:),xloc(2,:));
                Uds(3,:) = griddata(xloc_radial_grid(1,:),xloc_radial_grid(2,:),Uds_radial_grid(3,:),xloc(1,:),xloc(2,:));
                Uds(isnan(Uds))=0;
         
                dy = xystats(2,2)-xystats(2,1)
                gradient(Uds(1,:),
                Dx =  (Uds_plusx - Uds_minusx)/(2*dx);
                Dy =  (Uds_plusy - Uds_minusy)/(2*dy);
                
                Exx_ds = Dx(1,:)';
                Exy_ds = .5*(Dx(2,:)+Dy(1,:))';
                Eyy_ds = Dy(2,:)';
                
                %add elastic backslip contribution
                 [Uss,Dds] = backslip_elastic3D(m_ds,xloc);
             
                   Exx_ds = Exx_ds + Dds(1,:)';
                   Exy_ds = Exy_ds + .5*(Dds(2,:)+Dds(4,:))';
                   Eyy_ds = Eyy_ds + Dds(5,:)';
             
                   Exx_ds_elastic =  Dds(1,:)';
                   Exy_ds_elastic = .5*(Dds(2,:)+Dds(4,:))';
                   Eyy_ds_elastic = Dds(5,:)';
             
                
            else
                
                Exx_ds = zeros(size(xloc,2),1);
                Exy_ds = zeros(size(xloc,2),1);
                Eyy_ds = zeros(size(xloc,2),1);
                
                Exx_ds_elastic = Exx_ds;
                Exy_ds_elastic = Exy_ds;
                Eyy_ds_elastic = Eyy_ds;
                
            end
        
        
     end %if isnan(T(k))
     
  
     
    Gss.Exx(:,k) = Exx_ss;
    Gss.Exy(:,k) = Exy_ss;
    Gss.Eyy(:,k) = Eyy_ss;
  
    Gds.Exx(:,k) = Exx_ds;
    Gds.Exy(:,k) = Exy_ds;
    Gds.Eyy(:,k) = Eyy_ds;
  
    Gss_e.Exx(:,k) = Exx_ss;
    Gss_e.Exy(:,k) = Exy_ss;
    Gss_e.Eyy(:,k) = Eyy_ss;
  
    Gds_e.Exx(:,k) = Exx_ds_elastic;
    Gds_e.Exy(:,k) = Exy_ds_elastic;
    Gds_e.Eyy(:,k) = Eyy_ds_elastic;
    

     
end %for k=1:size(centers,1)


Gss.Exx = sum(Gss.Exx,2);
Gss.Exy = sum(Gss.Exy,2);
Gss.Eyy = sum(Gss.Eyy,2);

Gds.Exx = sum(Gds.Exx,2);
Gds.Exy = sum(Gds.Exy,2);
Gds.Eyy = sum(Gds.Eyy,2);

Gss_e.Exx = sum(Gss_e.Exx,2);
Gss_e.Exy = sum(Gss_e.Exy,2);
Gss_e.Eyy = sum(Gss_e.Eyy,2);

Gds_e.Exx = sum(Gds_e.Exx,2);
Gds_e.Exy = sum(Gds_e.Exy,2);
Gds_e.Eyy = sum(Gds_e.Eyy,2);
