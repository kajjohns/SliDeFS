function [GExx,GExy,GEyy]=make_StrainG_cycles(pm,rake_patches,xystats,H1,H2,tR1,tR2,Ts,teqs,Nterms)


%% specify regular grid of points to compute viscoelastic part of solution
xs = linspace(min(xystats(:,1))-100,max(xystats(:,1))+100,150);
ys = linspace(min(xystats(:,2)),max(xystats(:,2)),100);
[Xs,Ys] = meshgrid(xs,ys);
xy_visco = [Xs(:) Ys(:)];

for k=1:size(pm,1)
        
  [Exx,Exy,Eyy] = Get_interseismic_strainrates_cycle(pm(k,:),rake_patches(k),...
        xystats,xs,ys,xy_visco,H1,H2,tR1,tR2,Ts(k),teqs(k),Nterms);


        GExx(:,k) = Exx;
        GExy(:,k) = Exy;
        GEyy(:,k) = Eyy;

    disp(['completed ' num2str(k/size(pm,1)*100) ' % of viscoelastic calculations'])
      
end

