%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%this script plots moment source directions and uncertainty wedges
%Note: run sample_posterior first
%
%

vecscale = 5;  %scale for plotting vectors

%sizse of scale bar for legend, Nm/yr
legend_scale = 2*10^16; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%number of slip patches
numslip = size(pms,1);
all_moments_Ms = all_Ms(1+numslip:end,:); %extract moment rates (exclude slip sources)

%extract moment components
m11s = all_moments_Ms(1:end/3,:);
m12s = all_moments_Ms(1+end/3:2*end/3,:);
m22s = all_moments_Ms(1+2*end/3:end,:);


%% compute principal directions
minVecs = zeros(2,size(m11s,1),size(m11s,2));
maxVecs = zeros(2,size(m11s,1),size(m11s,2));
minvals = zeros(size(m11s,1),size(m11s,2));
maxvals = zeros(size(m11s,1),size(m11s,2));


for j=1:size(m11s,2)

        for k=1:size(m11s(:,j),1)

            
            E = [m11s(k,j) m12s(k,j); m12s(k,j) m22s(k,j)];
            [vec,val] = eig(E);

            minVecs(:,k,j) = vec(:,1);
            maxVecs(:,k,j) = vec(:,2);
            minvals(k,j) = val(1,1);
            maxvals(k,j) = val(2,2);
            
        end
        
end



mean_minvals = mean(minvals,2);
mean_maxvals = mean(maxvals,2);


f=figure;
hold on

skip=1;

for k=1:length(nodes)
    
    
        
    %plot the larger magnitude principal direction
    
    if abs(mean_minvals(k))>abs(mean_maxvals(k))
        
        bigvals = mean_minvals(k);
        smallvals = mean_maxvals(k);
        
        bigVecs = squeeze(minVecs(:,k,:));
        smallVecs = squeeze(maxVecs(:,k,:));
        

                
    else
        
        bigvals = mean_maxvals(k);
        smallvals = mean_minvals(k);        
    
        bigVecs = squeeze(maxVecs(:,k,:));
        smallVecs = squeeze(minVecs(:,k,:));

    end
    
    mean_smallVecs = mean(smallVecs,2);
    angles = squeeze(atan(bigVecs(2,:)./bigVecs(1,:))*180/pi);
    
    
    %determin if bimodal. If so, add 180 to negative values
    if ~isnan(sum(angles))
    [dip, p_value, xlow, xup] = HartigansDipSignifTest(angles, 50);
    if p_value<0.5 & mean(abs(angles))>45
        angles(angles<0) = 180 + angles(angles<0);
    end
    end
    
    mean_angles(k) = mean(angles);
    
    std_angles = std(angles);
    
    if mod(k,skip)==0
        
   
        
        a1 = pi/180*(mean(angles) + 2*std_angles);
        a2 = pi/180*(mean(angles) - 2*std_angles);
        t = linspace(a1,a2,128);



        x0 = [0 vecscale*cos(t) 0]*bigvals+nodes_llh(k,1);
        z0 = [0 vecscale*sin(t) 0]*bigvals+nodes_llh(k,2);    
        if bigvals<0
             patch( x0, z0, 'r','FaceAlpha',.25 ,'EdgeColor','none');
        else
             patch( x0, z0, 'b','FaceAlpha',.25 ,'EdgeColor','none');
        end    

        x0 = [0 -vecscale*cos(t) 0]*bigvals+nodes_llh(k,1);
        z0 = [0 -vecscale*sin(t) 0]*bigvals+nodes_llh(k,2);

        if bigvals<0
           patch( x0, z0, 'r','FaceAlpha',.25 ,'EdgeColor','none');
        else
            patch( x0, z0, 'b','FaceAlpha',.25 ,'EdgeColor','none');
        end



        %plot smaller principal direction as a  line
        vx = [nodes_llh(k,1)-mean_smallVecs(1)*vecscale*smallvals nodes_llh(k,1)+mean_smallVecs(1)*vecscale*smallvals];
        vy = [nodes_llh(k,2)-mean_smallVecs(2)*vecscale*smallvals nodes_llh(k,2)+mean_smallVecs(2)*vecscale*smallvals];


        if smallvals<0
            plot(vx,vy,'r')
        else
            plot(vx,vy,'b')
        end

    end
   
    
end




plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('Principal Directions, Moment Sources, 95% wedges') 
set(gca,'fontsize',15)



%plot scale
%moments are unscaled:  to scale moments, multiply by conv where conv is
%computed in lines 379-381 of invert_strainrate_for_backslip.m

plot_legend = legend_scale/conv;  %size of legend scale bar

X = xlim;
Y = ylim;
dY = 0.1*(Y(2)-Y(1));
scale_pos = [X(1)+0.1*(X(2)-X(1)) Y(1)-dY];
vx = [scale_pos(1)-vecscale*plot_legend scale_pos(1)+vecscale*plot_legend];
vy = [scale_pos(2) scale_pos(2)];
plot(vx,vy,'r','linewidth',1)
plot(vx,vy-dY/2,'b','linewidth',1)
text(scale_pos(1)-vecscale*plot_legend,scale_pos(2)+dY/2,['Scale, ' num2str(legend_scale) ' Nm/yr'])
text(scale_pos(1)+2*vecscale*plot_legend,scale_pos(2),' Contractional')
text(scale_pos(1)+2*vecscale*plot_legend,scale_pos(2)-dY/2,' Extensional')
set(gca,'fontsize',15)
axis equal; axis tight
X = xlim; Y = ylim;
dX = 0.05*(Y(2)-Y(1)); dY = 0.05*(Y(2)-Y(1));
xlim([X(1)-dX X(2)])
ylim([Y(1)-dY Y(2)])

f.Position = [154    68   735   758];

