function [GExx,GExy,GEyy] = make_strainG(pm,rake,xystats)

%convert from km to m
xystats = 1000*xystats;
pm(:,1:3) = 1000*pm(:,1:3);
pm(:,6:7) = 1000*pm(:,6:7);

xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

GExx = zeros(size(xystats,1),npatches);
GExy = zeros(size(xystats,1),npatches);
GEyy = zeros(size(xystats,1),npatches);

pm(:,3) = pm(:,3)+1; %small shift to avoid physically impossible error

for k=1:npatches
    
    ss= cos(rake(k)*pi/180);
    ds = sin(rake(k)*pi/180);   
    
    m1=[pm(k,:) ss ds 0]';  %positive slip is R and LL sense
   
    [U1,D,S,flag]=disloc3d(m1,xloc,1,.25);


      
   
    Exx = D(1,:)';
    Exy = .5*(D(2,:)+D(4,:))';
    Eyy = D(5,:)';

   
    GExx(:,k) = Exx;
    GExy(:,k) = Exy; 
    GEyy(:,k) = Eyy;

     disp(['completed ' num2str(k/npatches*100) ' percent deep patches'])
     
end

