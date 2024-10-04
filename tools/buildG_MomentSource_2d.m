function [Gexx,Gexy,Geyy] = buildG_MomentSource_2d(xy_obs,tri,nodes,tri_centroids,nu)

%convert from plane strain to plane stress by transforming elastic
%properties
nu = nu/(1+nu);  %Note: no need to transform Lame constants because solution only depends on nu
%GFs depend only on ratios of Lame constants
mu = 1; 
lam = 2*nu/(1-2*nu);


for k=1:size(tri,1)

    nd = nodes(tri(k,:),:);

    [u1,u2,e11,e12,e22] = triangular_moment_source(nd,tri_centroids(:,1:2),nu);

    %NOTE:  the GFs are computed for m12+m21, so need to add together
    e11.m12 = e11.m12 + e11.m21;
    e12.m12 = e12.m12 + e12.m21;
    e22.m12 = e22.m12 + e22.m21;
    

    %interpolate from triangle centroids to observation coordinates
    Gexx.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m11,xy_obs(:,1),xy_obs(:,2)); 
    i = isnan(Gexx.m11(:,k));

    Gexx.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m11,xy_obs(:,1),xy_obs(:,2)); Gexx.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m11,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m11,xy_obs(:,1),xy_obs(:,2)); Gexy.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m11,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m11(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m11,xy_obs(:,1),xy_obs(:,2)); Geyy.m11(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m11,xy_obs(i,1),xy_obs(i,2),'nearest');

    Gexx.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m12,xy_obs(:,1),xy_obs(:,2)); Gexx.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m12,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m12,xy_obs(:,1),xy_obs(:,2)); Gexy.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m12,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m12(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m12,xy_obs(:,1),xy_obs(:,2)); Geyy.m12(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m12,xy_obs(i,1),xy_obs(i,2),'nearest');

    Gexx.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m22,xy_obs(:,1),xy_obs(:,2)); Gexx.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e11.m22,xy_obs(i,1),xy_obs(i,2),'nearest');
    Gexy.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m22,xy_obs(:,1),xy_obs(:,2)); Gexy.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e12.m22,xy_obs(i,1),xy_obs(i,2),'nearest');
    Geyy.m22(:,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m22,xy_obs(:,1),xy_obs(:,2)); Geyy.m22(i,k) = griddata(tri_centroids(:,1),tri_centroids(:,2),e22.m22,xy_obs(i,1),xy_obs(i,2),'nearest');


     disp(['Completed ' num2str(k/size(tri,1)*100) '% of moment source Greens function calculations']) 
     
end
