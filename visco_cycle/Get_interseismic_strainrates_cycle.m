function [Exx,Exy,Eyy] = Get_interseismic_strainrates_cycle(pm,rake,xystats,xs,ys,xy_visco,...
                      H1,H2,tR1,tR2,T,timeEq,Nterms)

    

xloc = [xystats';zeros(1,size(xystats,1))];
xloc_visco_grid = [xy_visco';zeros(1,size(xy_visco,1))];
          


    

if isnan(T)  %no viscoelastic calculation


        Exx = zeros(size(xloc,2),1);
        Exy = zeros(size(xloc,2),1);
        Eyy = zeros(size(xloc,2),1);


else  %compute viscoelastic and elastic parts


        m = [pm cos(rake*pi/180) sin(rake*pi/180) 0];

        U=backslip_cycle3D_layered_noelastic(m,xloc_visco_grid,H1,H2,1,1,timeEq,tR1,tR2,T,Nterms);


        Ue = reshape(U(1,:),length(ys),length(xs));
        Un = reshape(U(2,:),length(ys),length(xs));


        [Due_Dx,Due_Dy] = gradient(Ue,xs,ys);
        [Dun_Dx,Dun_Dy] = gradient(Un,xs,ys);

        Exx_ss = Due_Dx(:);
        Exy_ss = .5*(Due_Dy(:)+Dun_Dx(:));
        Eyy_ss = Dun_Dy(:);


        %interpolate to desired grid
        Exx = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),Exx_ss,xloc(1,:),xloc(2,:))';
        Exy = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),Exy_ss,xloc(1,:),xloc(2,:))';
        Eyy = griddata(xloc_visco_grid(1,:),xloc_visco_grid(2,:),Eyy_ss,xloc(1,:),xloc(2,:))';
        Exx(isnan(Exx))=0;
        Exy(isnan(Exy))=0;
        Eyy(isnan(Eyy))=0;


end


               
                
         