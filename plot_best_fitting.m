

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
plot_earthquake_timing = true;

%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_max_shear
        %plot maximum shear strain rate
        f=figure;
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
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        f=figure;
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_bs,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_bs,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('model mean maximum shear strain rate (backslip)')
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)



        f=figure;
        if regular_grid 
            X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(max_shear_offfault,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,max_shear_offfault,'fill')
        end
        colormap(jet)
        colorbar
        set(gca,'colorscale','log')
        caxis([1e-3 1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('observed  - model max shear strain rate ')
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

end



if plot_dilatation
    
        load cmap
        f=figure;
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
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)


        f=figure;
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
        title('model dilatation rate')
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)



        f=figure;
        if regular_grid 
           X = [obs_llh(:,1)-dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)+dl/2 obs_llh(:,1)-dl/2];
            Y = [obs_llh(:,2)-dl/2 obs_llh(:,2)-dl/2 obs_llh(:,2)+dl/2 obs_llh(:,2)+dl/2];
            patch(X',Y',repmat(dilatation-dilatation_bs,1,4)','edgecolor','none');
        else
            scatter(obs_llh(:,1),obs_llh(:,2),80,dilatation-dilatation_bs,'fill')
        end
        colormap(cmap)
        colorbar
        caxis([-1e-1 1e-1])
        hold on
        plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
        title('observed - model dilatation rate')
        f.Position = [154    68   735   758];
        set(gca,'fontsize',15)

end

if plot_slip_deficit
    
        %plot slip deficit rate
        figure
        plotpatchslip3D_noedge(pm_top,bs_top*1000,1)
        hold on
        plotpatchslip3D_noedge(pm,bs_bot*1000,1)
        title('slip deficit rates, mm/yr')
        axis tight
        view(-35,35)
        set(gca,'fontsize',15)

        %plot slip deficit rate
        figure
        plotpatchslip3D_noedge(pm_top,bs_top*1000,1)
        hold on
        plotpatchslip3D_noedge(pm,bs_bot*1000,1)
        title('slip deficit rates, mm/yr')
        axis tight
        view(-35,35)
        set(gca,'colorscale','log')
        set(gca,'fontsize',15)

end


if plot_rake
    
        %plot rake
        figure
        plotpatchslip3D_noedge(pm_top,rake_top,1)
        hold on
        plotpatchslip3D_noedge(pm,rake_bot,1)
        title('rake, degrees')
        axis tight
        view(-35,35)
        set(gca,'fontsize',15)
        colormap(pmkmp(128,'IsoAZ'))
        caxis([-180 180])


        %plot preferred (prior) rake
        figure
        plotpatchslip3D_noedge(pm_top,rake_top_patches,1)
        hold on
        plotpatchslip3D_noedge(pm,rake_patches,1)
        title('pior (preferred) rake, degrees')
        axis tight
        view(-35,35)
        set(gca,'fontsize',15)
        colormap(pmkmp(128,'IsoAZ'))
        caxis([-180 180])

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



        figure
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
        c.FontSize=15
        title('depth-averaged slip deficit rates (mm/yr)','fontsize',15)    

        figure
        hold on
        for k=1:size(pm_ends,1)
            cline([pm_ends(k,1) pm_ends(k,3)],[pm_ends(k,2) pm_ends(k,4)],[depth_averaged_s(k) depth_averaged_s(k)]*1000)  %convert to mm/yr 
        end
        colormap(jet)
        caxis([0 40])
        c=colorbar;
        axis equal
        axis tight
        c.FontSize=15
        title('depth-averaged slip deficit rates (mm/yr)','fontsize',15)    


end

if plot_earthquake_timing

        figure
        plotpatchslip3D_noedge(pm_top,Ts_top,1)
        hold on
        plotpatchslip3D_noedge(pm,Ts,1)
        title('earthquake recurrence time, yr')
        axis tight
        view(-20,60)



        figure
        plotpatchslip3D_noedge(pm_top,teqs_top,1)
        hold on
        plotpatchslip3D_noedge(pm,teqs,1)
        title('time since last earthquake, yr')
        axis tight
        view(-20,60)
        
end
