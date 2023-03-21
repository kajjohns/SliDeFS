
%this script plots slip rates along faults
%Note: run sample_posterior first

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SECTION

%list of fault names to plot (note, only part of the name needs to be in
%the list to include the fault 
plot_names = {'San Andreas','Calaveras','Hosgri'};


%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for k=1:length(FaultInfo)
FaultID(k) = FaultInfo{k}.FaultID;
FaultNames{k} = FaultInfo{k}.FaultName;
end

%need to assign fault names to segments 
for k=1:size(SegEnds,1)   
    ind = segid(k) == FaultID;    
    fault_name{k} =  FaultNames{ind};  
end

%compute depth-averaged slip rate samples
top_row_samples_s = all_Ms(1:end/nve,:);
bottom = all_Ms(1+end/nve:end,:);
sum_bottom_rows = zeros(size(top_row_samples_s));
for k=1:nve-1
    sum_bottom_rows = sum_bottom_rows + bottom(k:nve-1:end,:);
end
depth_averaged_samples_s = 1/nve*(top_row_samples_s + sum_bottom_rows);

if variable_rake
    bs = sqrt(mhat_rake.^2 + mhat_perp.^2);
else
    bs = mhat_rake;
end
%plot depth-averaged slip deficit rate
top_row_s = bs(1:end/nve);
bottom = bs(1+end/nve:end);
for k=1:nve-1
    bottom_rows_s(:,k) = bottom(k:nve-1:end,:);
end
depth_averaged_s = 1/nve*(top_row_s + sum(bottom_rows_s,2));
        

if variable_rake
    depth_averaged_samples_s = sqrt(depth_averaged_samples_s(1:end/2,:).^2 + depth_averaged_samples_s(1+end/2:end,:).^2);
end

%get 95% range of depth-averaged slip deficit
s = depth_averaged_s*1000;  %convert to mm/yr
sorted = sort(depth_averaged_samples_s*1000,2); %convert to mm/yr
discard = round(.025*size(sorted,2));
sorted(:,1:discard) = [];
sorted(:,end-discard+1:end) = [];

upper_95 = sorted(:,end)-s;
lower_95 = s-sorted(:,1);

%because of poor sampling lower 95 can slip below zero -- fix that
lower_95(lower_95<0)=0;



for loop=1:length(plot_names)

    name = plot_names{loop};
    i=strmatch(name,fault_name);
    j=ismember(pm_top_seg_num,i);
    srate = s(j);

    upper_bound = UB(j);
    L95 = lower_95(j);
    U95 = upper_95(j);
    upper_bound_cfm = ub_top_patches(j);
   % pref = mean_s_prior(j);



    pm_x = pm_top(j,6);
    pm_y = pm_top(j,7);
    [y,i] = max(pm_y);
    first_x = pm_x(i);
    first_y = pm_y(i);
    
    distance_along_fault = sqrt((pm_x-first_x).^2 + (pm_y-first_y).^2);
    [y,i] = sort(distance_along_fault);


    figure
    subplot(121)
    phand(1) = errorbar(distance_along_fault(i),srate(i),L95(i),U95(i),'k','linewidth',2);
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    xlabel('distance from north end of fault, km')
    ylabel('depth-averaged slip deficit rate, mm/yr')
    title([name ', depth-averaged slip deficit'],'fontsize',15)
    
    hold on
    plot(distance_along_fault(i),upper_bound_cfm(i),'b','linewidth',3)
    plot(distance_along_fault(i),upper_bound(i),'r','linewidth',3)
   % plot(distance_along_fault(i),pref(i),'r:','linewidth',3)
    
    grid on
    legend('estimated','upper bound','lower bound')
    set(gca,'fontsize',15)
    
    subplot(122)
    hold on
    plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k')
    plot(pm_x(i),pm_y(i),'r','linewidth',3)
    set(gca,'fontsize',15)
    
    set(gcf, 'Position', get(0, 'Screensize'));

end


