

%pm_ID = GeoBounds(pm_top_seg_num,1);
pm_ID = pm_top_seg_num;

%build Gcreep -- piecewise linear smoothing along faults
faultnums = unique(floor(pm_ID));
cnt = 0;
GExx_top_fine = zeros(size(G1Exx_top,1),size(pm_top,1));
GExy_top_fine = zeros(size(G1Exx_top,1),size(pm_top,1));
GEyy_top_fine = zeros(size(G1Exx_top,1),size(pm_top,1));

for k=1:length(faultnums)

   n = sum(floor(pm_ID)==faultnums(k));
   
   if n==1
   
        
        GExx_top_fine(:,cnt+1) = G1Exx_top(:,cnt+1);  
        GExy_top_fine(:,cnt+1) = G1Exy_top(:,cnt+1);  
        GEyy_top_fine(:,cnt+1) = G1Eyy_top(:,cnt+1);  
       
   elseif n==2
       
       GExx_top_fine(:,cnt+1) = G1Exx_top(:,cnt+1)-G2Exx_top(:,cnt+1);  
       GExx_top_fine(:,cnt+2) = G2Exx_top(:,cnt+1)+G1Exx_top(:,cnt+2);  
       
       GExy_top_fine(:,cnt+1) = G1Exy_top(:,cnt+1)-G2Exy_top(:,cnt+1);  
       GExy_top_fine(:,cnt+2) = G2Exy_top(:,cnt+1)+G1Exy_top(:,cnt+2);  
       
       GEyy_top_fine(:,cnt+1) = G1Eyy_top(:,cnt+1)-G2Eyy_top(:,cnt+1);  
       GEyy_top_fine(:,cnt+2) = G2Eyy_top(:,cnt+1)+G1Eyy_top(:,cnt+2);  

   else
 
        GExx_top_fine(:,cnt+1) = G1Exx_top(:,cnt+1)-G2Exx_top(:,cnt+1);  
        GExy_top_fine(:,cnt+1) = G1Exy_top(:,cnt+1)-G2Exy_top(:,cnt+1);  
        GEyy_top_fine(:,cnt+1) = G1Eyy_top(:,cnt+1)-G2Eyy_top(:,cnt+1);  
 
       for j=cnt+2:cnt+(n-1)
             GExx_top_fine(:,j) = G2Exx_top(:,j-1)+G1Exx_top(:,j)-G2Exx_top(:,j);  
            GExy_top_fine(:,j) = G2Exy_top(:,j-1)+G1Exy_top(:,j)-G2Exy_top(:,j);  
            GEyy_top_fine(:,j) = G2Eyy_top(:,j-1)+G1Eyy_top(:,j)-G2Eyy_top(:,j);  
            
       end     
       
         GExx_top_fine(:,cnt+n) = G2Exx_top(:,cnt+n-1)+G1Exx_top(:,cnt+n);  
        GExy_top_fine(:,cnt+n) = G2Exy_top(:,cnt+n-1)+G1Exy_top(:,cnt+n);  
        GEyy_top_fine(:,cnt+n) = G2Eyy_top(:,cnt+n-1)+G1Eyy_top(:,cnt+n);  
 
             
   end
   
   

   cnt = cnt+n;
   
end

