function [G1Exx,G2Exx,G1Exy,G2Exy,G1Eyy,G2Eyy] = make_strainG_piecewise(pm,rake,xystats,refine)

%convert from km to m
xystats = 1000*xystats;
pm(:,1:3) = 1000*pm(:,1:3);
pm(:,6:7) = 1000*pm(:,6:7);
refine = 1000*refine;


xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

G1Exx = zeros(size(xystats,1),npatches);
G2Exx = zeros(size(xystats,1),npatches);
G1Exy = zeros(size(xystats,1),npatches);
G2Exy = zeros(size(xystats,1),npatches);
G1Eyy = zeros(size(xystats,1),npatches);
G2Eyy = zeros(size(xystats,1),npatches);


for k=1:npatches
    
    %divide into small segments 
    nhe = ceil(pm(k,1)/refine);
    pf = patchfault(pm(k,:),nhe,1);
    pf(:,3) = pf(:,3)+1; %small shift to avoid physically impossible error
    
    Exx1 = zeros(size(xystats,1),1);
    Exx2 = zeros(size(xystats,1),1);
    Exy1 = zeros(size(xystats,1),1);
    Exy2 = zeros(size(xystats,1),1);
    Eyy1 = zeros(size(xystats,1),1);
    Eyy2 = zeros(size(xystats,1),1);

    ss= cos(rake(k)*pi/180);
    ds = sin(rake(k)*pi/180);

    
    for j=1:nhe
        
        m1=[pf(j,:) ss ds 0]';
   
        [U1,D,S,flag]=disloc3d(m1,xloc,1,.25);
       
        
        Exx = D(1,:)';
        Exy = .5*(D(2,:)+D(4,:))';
        Eyy = D(5,:)';

        Exx1 = Exx1 + Exx;
        Exy1 = Exy1 + Exy;
        Eyy1 = Eyy1 + Eyy;
        

        Exx2 = Exx2 + (j-1)/nhe*Exx;
        Exy2 = Exy2 + (j-1)/nhe*Exy;
        Eyy2 = Eyy2 + (j-1)/nhe*Eyy;

        
    end
    
  

    G1Exx(:,k) = Exx1;
    G1Exy(:,k) = Exy1; 
    G1Eyy(:,k) = Eyy1;

    G2Exx(:,k) = Exx2;
    G2Exy(:,k) = Exy2; 
    G2Eyy(:,k) = Eyy2;

    
    disp(['Completed ' num2str(k/npatches*100) ' % of shallow patch Greens Functions'])
end

