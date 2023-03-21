
list_id = unique(pm_top_seg_id);
N = size(pm_top,1);

%horizontal and vertical differences
pm_smooth_hz = zeros(N+size(pm,1),N+size(pm,1));
pm_smooth_vert = zeros(N+size(pm,1),N+size(pm,1));

for k=1:length(list_id)
   
    i = pm_top_seg_id == list_id(k);
    
    ps = find(i);
    for j=1:length(ps)
        
        
        %top row
       
        
        if j<length(ps)
             pm_smooth_hz(ps(j),ps(j)) = 1;
             pm_smooth_hz(ps(j),ps(j)+1) = -1; %smooth patch to left
        end
        
        pm_smooth_vert(ps(j),ps(j)) = 1;
        pm_smooth_vert(ps(j),N+(nve-1)*(ps(j)-1)+1) = -1;  %smooth patch below
        
        
        %next rows
        for m = 1:nve-2
            
            
             if j<length(ps)
                 pm_smooth_hz(N+(nve-1)*(ps(j)-1)+m,N+(nve-1)*(ps(j)-1)+m) = 1; 
                 pm_smooth_hz(N+(nve-1)*(ps(j)-1)+m,N+(nve-1)*(ps(j)-1)+m+(nve-1)) = -1;
             end
             
            pm_smooth_vert(N+(nve-1)*(ps(j)-1)+m,N+(nve-1)*(ps(j)-1)+m) = 1; 
            pm_smooth_vert(N+(nve-1)*(ps(j)-1)+m,N+(nve-1)*(ps(j)-1)+m+1) = -1; 
            
        end
        
    
    end
     
    
    
end


%remove zero rows
pm_smooth_hz(sum(abs(pm_smooth_hz),2)==0,:)=[];
pm_smooth_vert(sum(abs(pm_smooth_vert),2)==0,:)=[];

pm_smooth = [pm_smooth_hz;pm_smooth_vert];
