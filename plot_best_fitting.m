

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SECTION

%select plots to make (true/false)

%map-view plots of observed and model strain rates
plot_max_shear = true;
plot_dilatation = true;

%plot slip deficit on slip patches
plot_slip_deficit = true;

%plot rake on slip patches
plot_rake = true;

%plot depth-averaged slip rates along fault traces
plot_depth_averaged = true;

%plot recurrence times and times since last earthquake
plot_earthquake_timing = false;

%plot distributed moment sources
plot_moment = true;


%scale the size of strain rate axes for plotting
vecscale = .5;

%scale the size of moment rate axes for plotting
vecscale_mom = 1;


%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_max_shear
        %plot maximum shear strain rate
        f=figure;
        subplot(211)
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('observed mean maximum shear strain rate')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        subplot(212)
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_bs,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_total,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('model mean maximum shear strain rate (total)')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        

        f=figure;
        subplot(211)
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_offfault,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_bs,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('max shear strain rate due to backslip')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        if include_moment 

        subplot(212)
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_offfault,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_bf,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('max shear strain rate due to moment sources')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

        end %if include_moment 

        f=figure;
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_offfault,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_residual,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('residual max shear strain rate ')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

end



if plot_dilatation
    
        load cmap
        f=figure;
        subplot(211)          
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('observed dilatation rate')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

        subplot(212)          
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation_total,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation_total,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('model dilatation rate (total)')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        f=figure;
        subplot(211)          
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation_bs,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation_bs,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('dilatation rate due to backslip')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

        subplot(212)          
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation_bf,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation_bf,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('dilatation rate due to moment sources')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        f=figure;
        if regular_grid 
           X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation_residual,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation_residual,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('residual dilatation rate')
        axis equal; axis tight
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

end

if plot_slip_deficit
    
        %plot slip deficit rate
        f=figure;
        plotpatchslip3D_noedge(pm_top,bs_top*1000,1)
        hold on
        plotpatchslip3D_noedge(pm,bs_bot*1000,1)
        title('slip deficit rates, mm/yr')
        axis equal; axis tight
        view(-35,35)
        set(gca,'fontsize',15)
        f.Position = [154    68   735   758];

        %plot slip deficit rate
        f=figure;
        plotpatchslip3D_noedge(pm_top,bs_top*1000,1)
        hold on
        plotpatchslip3D_noedge(pm,bs_bot*1000,1)
        title('slip deficit rates, mm/yr')
        axis equal; axis tight
        view(-35,35)
        set(gca,'colorscale','log')
        set(gca,'fontsize',15)
        f.Position = [154    68   735   758];

end


if plot_rake
    
        %plot rake
        f=figure;
        plotpatchslip3D_noedge(pm_top,rake_top,1)
        hold on
        plotpatchslip3D_noedge(pm,rake_bot,1)
        title('rake, degrees')
        axis equal; axis tight
        view(-35,35)
        set(gca,'fontsize',15)
        colormap(pmkmp(128,'IsoAZ'))
        caxis([-180 180])
        f.Position = [154    68   735   758];

        %plot preferred (prior) rake
        f=figure;
        plotpatchslip3D_noedge(pm_top,rake_top_patches,1)
        hold on
        plotpatchslip3D_noedge(pm,rake_patches,1)
        title('pior (preferred) rake, degrees')
        axis equal; axis tight
        view(-35,35)
        set(gca,'fontsize',15)
        colormap(pmkmp(128,'IsoAZ'))
        caxis([-180 180])
        f.Position = [154    68   735   758];

end


if plot_depth_averaged
        
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

        a = 90-pm_top(:,5);
        pm_ends = [pm_top(:,6)+pm_top(:,1).*cos(a*pi/180)/2  pm_top(:,7)+pm_top(:,1).*sin(a*pi/180)/2  pm_top(:,6)-pm_top(:,1).*cos(a*pi/180)/2  pm_top(:,7)-pm_top(:,1).*sin(a*pi/180)/2];



        f=figure;
        hold on
        for k=1:size(pm_ends,1)
            cline([pm_ends(k,1) pm_ends(k,3)],[pm_ends(k,2) pm_ends(k,4)],[depth_averaged_s(k) depth_averaged_s(k)]*1000)  %convert to mm/yr 
        end
        colormap(jet)
        set(gca,'ColorScale','log')
        caxis([0.1 40])
        c=colorbar;
        axis equal
        axis tight
        c.FontSize=15;
        f.Position = [154    68   735   758];
        title('depth-averaged slip deficit rates (mm/yr)','fontsize',15)    

        f=figure;
        hold on
        for k=1:size(pm_ends,1)
            cline([pm_ends(k,1) pm_ends(k,3)],[pm_ends(k,2) pm_ends(k,4)],[depth_averaged_s(k) depth_averaged_s(k)]*1000)  %convert to mm/yr 
        end
        colormap(jet)
        caxis([0 40])
        c=colorbar;
        axis equal
        axis tight
        c.FontSize=15;
        f.Position = [154    68   735   758];
        title('depth-averaged slip deficit rates (mm/yr)','fontsize',15)    


end

if plot_earthquake_timing

        f=figure;
        plotpatchslip3D_noedge(pm_top,Ts_top,1)
        hold on
        plotpatchslip3D_noedge(pm,Ts,1)
        title('earthquake recurrence time, yr')
        axis equal; axis tight
        view(-20,60)
        f.Position = [154    68   735   758];



        f=figure;
        plotpatchslip3D_noedge(pm_top,teqs_top,1)
        hold on
        plotpatchslip3D_noedge(pm,teqs,1)
        title('time since last earthquake, yr')
        axis equal; axis tight
        view(-20,60)
        f.Position = [154    68   735   758];
        
end



%plot principal directions

vecscale = vecscale*10^6;  %for plotting 

%total
clear minVecs maxVecs minvals maxvals
for j=1:size(Exx_total,1)

    E = [Exx_total(j) Exy_total(j); Exy_total(j) Eyy_total(j)];
    [vec,val] = eig(E);

    minVecs(:,j) = vec(:,1);
    maxVecs(:,j) = vec(:,2);
    minvals(j) = val(1,1);
    maxvals(j) = val(2,2);
          
    
end

style_total = -(minvals+maxvals)./(abs(minvals)+abs(maxvals));

f=figure;
hold on

%plot model directions


for k=1:length(Exx_total)
    
    
    %plot the larger magnitude principal direction
    
    if abs(minvals(k))>abs(maxvals(k))
        
        bigvals = minvals(k);
        smallvals = maxvals(k);
        
        bigVecs = minVecs(:,k);
        smallVecs = maxVecs(:,k);
        

                
    else
        
        bigvals = maxvals(k);
        smallvals = minvals(k);        
    
        bigVecs = maxVecs(:,k);
        smallVecs = minVecs(:,k);

    end
       
    
     angles(k) = atan(bigVecs(2,:)./bigVecs(1,:))*180/pi;
   
       

            %plot larger principal direction 
            vx = [obs_llh(k,1)-bigVecs(1)*bigvals*vecscale obs_llh(k,1)+bigVecs(1)*bigvals*vecscale];
            vy = [obs_llh(k,2)-bigVecs(2)*bigvals*vecscale obs_llh(k,2)+bigVecs(2)*bigvals*vecscale];


            if bigvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end


            %plot smaller principal direction 
            vx = [obs_llh(k,1)-smallVecs(1)*vecscale*smallvals obs_llh(k,1)+smallVecs(1)*vecscale*smallvals];
            vy = [obs_llh(k,2)-smallVecs(2)*vecscale*smallvals obs_llh(k,2)+smallVecs(2)*vecscale*smallvals];


            if smallvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end

        
         

 end
    

plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)

%plot scale
X = xlim;
Y = ylim;
dY = 0.1*(Y(2)-Y(1));
scale_pos = [X(1)+0.1*(X(2)-X(1)) Y(1)-dY];
scaleval=10^-7;
vx = [scale_pos(1)-vecscale*scaleval scale_pos(1)+vecscale*scaleval];
vy = [scale_pos(2) scale_pos(2)];
plot(vx,vy,'r','linewidth',1)
plot(vx,vy-dY/2,'b','linewidth',1)
text(scale_pos(1)-vecscale*scaleval,scale_pos(2)+dY/2,['Scale, ' num2str(scaleval) ' strain/yr'])
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2),' Contractional')
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2)-dY/2,' Extensional')
set(gca,'fontsize',15)
axis equal
axis tight
X = xlim; Y = ylim;
dX = 0.05*(Y(2)-Y(1)); dY = 0.05*(Y(2)-Y(1));
xlim([X(1)-dX X(2)])
ylim([Y(1)-dY Y(2)])
f.Position = [154    68   735   758];
title('Principal directions -- total strain rate')



%backslip
clear minVecs maxVecs minvals maxvals

for j=1:size(Exx_bs,1)

    E = [Exx_bs(j) Exy_bs(j); Exy_bs(j) Eyy_bs(j)];
    [vec,val] = eig(E);

    minVecs(:,j) = vec(:,1);
    maxVecs(:,j) = vec(:,2);
    minvals(j) = val(1,1);
    maxvals(j) = val(2,2);
          
        
end

style_bs = -(minvals+maxvals)./(abs(minvals)+abs(maxvals));



f=figure;

hold on


%plot model directions


for k=1:length(Exx_bs)
    
    
    %plot the larger magnitude principal direction
    
    if abs(minvals(k))>abs(maxvals(k))
        
        bigvals = minvals(k);
        smallvals = maxvals(k);
        
        bigVecs = minVecs(:,k);
        smallVecs = maxVecs(:,k);
        

                
    else
        
        bigvals = maxvals(k);
        smallvals = minvals(k);        
    
        bigVecs = maxVecs(:,k);
        smallVecs = minVecs(:,k);

    end
       
    
     angles(k) = atan(bigVecs(2,:)./bigVecs(1,:))*180/pi;
   
       

            %plot larger principal direction 
            vx = [obs_llh(k,1)-bigVecs(1)*bigvals*vecscale obs_llh(k,1)+bigVecs(1)*bigvals*vecscale];
            vy = [obs_llh(k,2)-bigVecs(2)*bigvals*vecscale obs_llh(k,2)+bigVecs(2)*bigvals*vecscale];


            if bigvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end


            %plot smaller principal direction 
            vx = [obs_llh(k,1)-smallVecs(1)*vecscale*smallvals obs_llh(k,1)+smallVecs(1)*vecscale*smallvals];
            vy = [obs_llh(k,2)-smallVecs(2)*vecscale*smallvals obs_llh(k,2)+smallVecs(2)*vecscale*smallvals];


            if smallvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end

        
         
     

 end
    

plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)

%plot scale
X = xlim;
Y = ylim;
dY = 0.1*(Y(2)-Y(1));
scale_pos = [X(1)+0.1*(X(2)-X(1)) Y(1)-dY];
scaleval=10^-7;
vx = [scale_pos(1)-vecscale*scaleval scale_pos(1)+vecscale*scaleval];
vy = [scale_pos(2) scale_pos(2)];
plot(vx,vy,'r','linewidth',1)
plot(vx,vy-dY/2,'b','linewidth',1)
text(scale_pos(1)-vecscale*scaleval,scale_pos(2)+dY/2,['Scale, ' num2str(scaleval) ' strain/yr'])
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2),' Contractional')
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2)-dY/2,' Extensional')
set(gca,'fontsize',15)
axis equal
axis tight
X = xlim; Y = ylim;
dX = 0.05*(Y(2)-Y(1)); dY = 0.05*(Y(2)-Y(1));
xlim([X(1)-dX X(2)])
ylim([Y(1)-dY Y(2)])

f.Position = [154    68   735   758];
title('Principal strain rates --  backslip contribution')




%moment sources
clear minVecs maxVecs minvals maxvals

for j=1:size(Exx_bf,1)

    E = [Exx_bf(j) Exy_bf(j); Exy_bf(j) Eyy_bf(j)];
    [vec,val] = eig(E);

    minVecs(:,j) = vec(:,1);
    maxVecs(:,j) = vec(:,2);
    minvals(j) = val(1,1);
    maxvals(j) = val(2,2);
          
        
end

style_bf = -(minvals+maxvals)./(abs(minvals)+abs(maxvals));


f=figure;

hold on


%plot model directions


for k=1:length(Exx_bf)
    
    
    %plot the larger magnitude principal direction
    
    if abs(minvals(k))>abs(maxvals(k))
        
        bigvals = minvals(k);
        smallvals = maxvals(k);
        
        bigVecs = minVecs(:,k);
        smallVecs = maxVecs(:,k);
        

                
    else
        
        bigvals = maxvals(k);
        smallvals = minvals(k);        
    
        bigVecs = maxVecs(:,k);
        smallVecs = minVecs(:,k);

    end
       
    
     angles(k) = atan(bigVecs(2,:)./bigVecs(1,:))*180/pi;
   
       

            %plot larger principal direction 
            vx = [obs_llh(k,1)-bigVecs(1)*bigvals*vecscale obs_llh(k,1)+bigVecs(1)*bigvals*vecscale];
            vy = [obs_llh(k,2)-bigVecs(2)*bigvals*vecscale obs_llh(k,2)+bigVecs(2)*bigvals*vecscale];


            if bigvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end


            %plot smaller principal direction 
            vx = [obs_llh(k,1)-smallVecs(1)*vecscale*smallvals obs_llh(k,1)+smallVecs(1)*vecscale*smallvals];
            vy = [obs_llh(k,2)-smallVecs(2)*vecscale*smallvals obs_llh(k,2)+smallVecs(2)*vecscale*smallvals];


            if smallvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end

        
         
         

 end
    

plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)


%plot scale
X = xlim;
Y = ylim;
dY = 0.1*(Y(2)-Y(1));
scale_pos = [X(1)+0.1*(X(2)-X(1)) Y(1)-dY];
scaleval=10^-7;
vx = [scale_pos(1)-vecscale*scaleval scale_pos(1)+vecscale*scaleval];
vy = [scale_pos(2) scale_pos(2)];
plot(vx,vy,'r','linewidth',1)
plot(vx,vy-dY/2,'b','linewidth',1)
text(scale_pos(1)-vecscale*scaleval,scale_pos(2)+dY/2,['Scale, ' num2str(scaleval) ' strain/yr'])
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2),' Contractional')
text(scale_pos(1)+2*vecscale*scaleval,scale_pos(2)-dY/2,' Extensional')
set(gca,'fontsize',15)
axis equal
axis tight
X = xlim; Y = ylim;
dX = 0.05*(Y(2)-Y(1)); dY = 0.05*(Y(2)-Y(1));
xlim([X(1)-dX X(2)])
ylim([Y(1)-dY Y(2)])


f.Position = [154    68   735   758];
set(gca,'fontsize',15)
title('Principal strain rates -- distributed moment source contribution')



f=figure;
scatter(obs_llh(:,1),obs_llh(:,2),30,style_total,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)
axis equal; axis tight
f.Position = [154    68   735   758];

title('Strain Rate Style -- total model,red=reverse, green=ss, blue=normal')

f=figure;
scatter(obs_llh(:,1),obs_llh(:,2),30,style_bs,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)
axis equal; axis tight
f.Position = [154    68   735   758];
title('Strain Rate Style -- backslip,red=reverse, green=ss, blue=normal')

f=figure;
scatter(obs_llh(:,1),obs_llh(:,2),30,style_bf,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)
axis equal; axis tight
f.Position = [154    68   735   758];
title('Strain Rate Style -- moment contribution,red=reverse, green=ss, blue=normal')


if plot_moment & include_moment

%plot moment source principal directions


clear minVecs maxVecs minvals maxvals

for j=1:size(m11,1)

    E = [m11(j) .5*m12_m21(j); .5*m12_m21(j) m22(j)];
    [vec,val] = eig(E);

    minVecs(:,j) = vec(:,1);
    maxVecs(:,j) = vec(:,2);
    minvals(j) = val(1,1);
    maxvals(j) = val(2,2);
          
        
end

style_mom = -(minvals+maxvals)./(abs(minvals)+abs(maxvals));



f=figure;

hold on


%plot moment source directions (normalized by area)
vecscale_mom = vecscale_mom*10^-14;

for k=1:length(m11)
    
    
    %plot the larger magnitude principal direction
    
    if abs(minvals(k))>abs(maxvals(k))
        
        bigvals = minvals(k)/tri_areas(k);
        smallvals = maxvals(k)/tri_areas(k);
        
        bigVecs = minVecs(:,k)/tri_areas(k);
        smallVecs = maxVecs(:,k)/tri_areas(k);
        

                
    else
        
        bigvals = maxvals(k)/tri_areas(k);
        smallvals = minvals(k)/tri_areas(k);        
    
        bigVecs = maxVecs(:,k)/tri_areas(k);
        smallVecs = minVecs(:,k)/tri_areas(k);

    end
       
    
     angles(k) = atan(bigVecs(2,:)./bigVecs(1,:))*180/pi;
   
       

            %plot larger principal direction 
            vx = [tri_centroids_llh(k,1)-bigVecs(1)*bigvals*vecscale_mom tri_centroids_llh(k,1)+bigVecs(1)*bigvals*vecscale_mom];
            vy = [tri_centroids_llh(k,2)-bigVecs(2)*bigvals*vecscale_mom tri_centroids_llh(k,2)+bigVecs(2)*bigvals*vecscale_mom];


            if bigvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end


            %plot smaller principal direction 
            vx = [tri_centroids_llh(k,1)-smallVecs(1)*vecscale_mom*smallvals tri_centroids_llh(k,1)+smallVecs(1)*vecscale_mom*smallvals];
            vy = [tri_centroids_llh(k,2)-smallVecs(2)*vecscale_mom*smallvals tri_centroids_llh(k,2)+smallVecs(2)*vecscale_mom*smallvals];


            if smallvals<0
                plot(vx,vy,'r','linewidth',1)
            else
                plot(vx,vy,'b','linewidth',1)
            end

        
         
         

 end
    

plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)

%plot scale
X = xlim;
Y = ylim;
dY = 0.1*(Y(2)-Y(1));
scale_pos = [X(1)+0.1*(X(2)-X(1)) Y(1)-dY];
scaleval=10^13;
vx = [scale_pos(1)-vecscale_mom*scaleval scale_pos(1)+vecscale_mom*scaleval];
vy = [scale_pos(2) scale_pos(2)];
plot(vx,vy,'r','linewidth',1)
plot(vx,vy-dY/2,'b','linewidth',1)
text(scale_pos(1)-vecscale_mom*scaleval,scale_pos(2)+dY/2,['Scale, ' num2str(scaleval,'%E') ' (Nm/yr)/km^2'])
text(scale_pos(1)+2*vecscale_mom*scaleval,scale_pos(2),' Contractional')
text(scale_pos(1)+2*vecscale_mom*scaleval,scale_pos(2)-dY/2,' Extensional')
set(gca,'fontsize',15)
axis equal; axis tight
X = xlim; Y = ylim;
dX = 0.05*(Y(2)-Y(1)); dY = 0.05*(Y(2)-Y(1));
xlim([X(1)-dX X(2)])
ylim([Y(1)-dY Y(2)])

f.Position = [154    68   735   758];
title('Moment Sources (moment per area) -- Principal directions ')


f=figure;
scatter(tri_centroids_llh(:,1),tri_centroids_llh(:,2),30,style_mom,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
set(gca,'fontsize',15)
axis equal
axis tight
f.Position = [154    68   735   758];
title('Moment Sources -- equiv. faulting style,red=reverse, green=ss, blue=normal ')


end
