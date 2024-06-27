function [Gexx,Gexy,Geyy] = buildG_MomentSource_2d(tri_centroids,xy_obs,nodes,nu)

%convert from plane strain to plane stress by transforming elastic
%properties
nu = nu/(1+nu);  %Note: no need to transform Lame constants because solution only depends on nu
%GFs depend only on ratios of Lame constants
mu = 1; 
lam = 2*nu/(1-2*nu);


for k=1:size(nodes,1)

    [df1dx,df2dx,df1dy,df2dy] = force_couple_planestrain(tri_centroids(:,1),tri_centroids(:,2),nodes(k,1),nodes(k,2),nu,mu);
    
    %force_coupl_planesrain returns stresses - need to to convert to strain
    %rate. 
    
    df1dx.exx = ((lam+2*mu)*df1dx.s11 - lam*df1dx.s22)/(2*mu*(lam+2*mu)+lam*(lam+2*mu)-lam^2);
    df1dx.exy = df1dx.s12/(2*mu);
    df1dx.eyy =  (df1dx.s22-lam*df1dx.exx)/(lam+2*mu);

    df2dx.exx = ((lam+2*mu)*df2dx.s11 - lam*df2dx.s22)/(2*mu*(lam+2*mu)+lam*(lam+2*mu)-lam^2);
    df2dx.exy = df2dx.s12/(2*mu);
    df2dx.eyy =  (df2dx.s22-lam*df2dx.exx)/(lam+2*mu);

    df1dy.exx = ((lam+2*mu)*df1dy.s11 - lam*df1dy.s22)/(2*mu*(lam+2*mu)+lam*(lam+2*mu)-lam^2);
    df1dy.exy = df1dy.s12/(2*mu);
    df1dy.eyy =  (df1dy.s22-lam*df1dy.exx)/(lam+2*mu);

    df2dy.exx = ((lam+2*mu)*df2dy.s11 - lam*df2dy.s22)/(2*mu*(lam+2*mu)+lam*(lam+2*mu)-lam^2);
    df2dy.exy = df2dy.s12/(2*mu);
    df2dy.eyy =  (df2dy.s22-lam*df2dy.exx)/(lam+2*mu);


    exx_m11 = df1dx.exx;
    exy_m11 = df1dx.exy;
    eyy_m11 = df1dx.eyy;

    exx_m12 = df1dy.exx + df2dx.exx;
    exy_m12 = df1dy.exy + df2dx.exy;
    eyy_m12 = df1dy.eyy + df2dx.eyy;

    
    exx_m22 = df2dy.exx;
    exy_m22 = df2dy.exy;
    eyy_m22 = df2dy.eyy;


   

    %interpolate from triangle centroids to observation coordinates
    Gexx.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m11,xy_obs(:,1),xy_obs(:,2)); 
    i = isnan(Gexx.m11(:,k));

    Gexx.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m11,xy_obs(:,1),xy_obs(:,2)); Gexx.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m11,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m11,xy_obs(:,1),xy_obs(:,2)); Gexy.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m11,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m11,xy_obs(:,1),xy_obs(:,2)); Geyy.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m11,xy_obs(i,1),xy_obs(i,2),'nearest');

    Gexx.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m12,xy_obs(:,1),xy_obs(:,2)); Gexx.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m12,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m12,xy_obs(:,1),xy_obs(:,2)); Gexy.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m12,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m12,xy_obs(:,1),xy_obs(:,2)); Geyy.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m12,xy_obs(i,1),xy_obs(i,2),'nearest');

    Gexx.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m22,xy_obs(:,1),xy_obs(:,2)); Gexx.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exx_m22,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m22,xy_obs(:,1),xy_obs(:,2)); Gexy.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),exy_m22,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m22,xy_obs(:,1),xy_obs(:,2)); Geyy.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),eyy_m22,xy_obs(i,1),xy_obs(i,2),'nearest');


     disp(['Completed ' num2str(k/size(nodes,1)*100) '% of moment source Greens function calculations']) 
     
end
