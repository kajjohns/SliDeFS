
pm_ID = SegID(pm_top_seg_num,1);

%build Gcreep -- piecewise linear smoothing along faults
faultnums = unique(floor(pm_ID));
cnt = 0;
GExx_top = zeros(size(G1Exx_top,1),size(pm_top,1));
GExy_top = zeros(size(G1Exx_top,1),size(pm_top,1));
GEyy_top = zeros(size(G1Exx_top,1),size(pm_top,1));

for k=1:length(faultnums)

   n = sum(floor(pm_ID)==faultnums(k));
   
   if n==1
   
        
        GExx_top(:,cnt+1) = G1Exx_top(:,cnt+1);  
        GExy_top(:,cnt+1) = G1Exy_top(:,cnt+1);  
        GEyy_top(:,cnt+1) = G1Eyy_top(:,cnt+1);  
       
   elseif n==2
       
       GExx_top(:,cnt+1) = G1Exx_top(:,cnt+1)-G2Exx_top(:,cnt+1);  
       GExx_top(:,cnt+2) = G2Exx_top(:,cnt+1)+G1Exx_top(:,cnt+2);  
       
       GExy_top(:,cnt+1) = G1Exy_top(:,cnt+1)-G2Exy_top(:,cnt+1);  
       GExy_top(:,cnt+2) = G2Exy_top(:,cnt+1)+G1Exy_top(:,cnt+2);  
       
       GEyy_top(:,cnt+1) = G1Eyy_top(:,cnt+1)-G2Eyy_top(:,cnt+1);  
       GEyy_top(:,cnt+2) = G2Eyy_top(:,cnt+1)+G1Eyy_top(:,cnt+2);  

   else
 
        GExx_top(:,cnt+1) = G1Exx_top(:,cnt+1)-G2Exx_top(:,cnt+1);  
        GExy_top(:,cnt+1) = G1Exy_top(:,cnt+1)-G2Exy_top(:,cnt+1);  
        GEyy_top(:,cnt+1) = G1Eyy_top(:,cnt+1)-G2Eyy_top(:,cnt+1);  
 
       for j=cnt+2:cnt+(n-1)
             GExx_top(:,j) = G2Exx_top(:,j-1)+G1Exx_top(:,j)-G2Exx_top(:,j);  
            GExy_top(:,j) = G2Exy_top(:,j-1)+G1Exy_top(:,j)-G2Exy_top(:,j);  
            GEyy_top(:,j) = G2Eyy_top(:,j-1)+G1Eyy_top(:,j)-G2Eyy_top(:,j);  
            
       end     
       
         GExx_top(:,cnt+n) = G2Exx_top(:,cnt+n-1)+G1Exx_top(:,cnt+n);  
        GExy_top(:,cnt+n) = G2Exy_top(:,cnt+n-1)+G1Exy_top(:,cnt+n);  
        GEyy_top(:,cnt+n) = G2Eyy_top(:,cnt+n-1)+G1Eyy_top(:,cnt+n);  
 
             
   end
   
   

   cnt = cnt+n;
   
end

