%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%No inputs required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%make realizations for computing uncertainties in nonlinear quantities
N = 1000;  %number of realizations
Exx_r = repmat(Exx_mean,1,N) +  repmat(Exx_std,1,N).*randn(size(Exx_mean,1),N);
Exy_r = repmat(Exy_mean,1,N) +  repmat(Exy_std,1,N).*randn(size(Exy_mean,1),N);
Eyy_r = repmat(Eyy_mean,1,N) +  repmat(Eyy_std,1,N).*randn(size(Eyy_mean,1),N);


max_shear_realizations = sqrt((Exx_r-Eyy_r).^2 + Exy_r.^2);
dilatation_realizations = Exx_r+Eyy_r;

std_maxshear = std(max_shear_realizations,[],2);
std_dilatation = std(dilatation_realizations,[],2);






%% compute principal directions
minVecs = zeros(2,size(Exx_r,1),size(Exx_r,2));
maxVecs = zeros(2,size(Exx_r,1),size(Exx_r,2));
minvals = zeros(size(Exx_r,1),size(Exx_r,2));
maxvals = zeros(size(Exx_r,1),size(Exx_r,2));


for j=1:size(Exx_r,2)

        for k=1:size(Exx_r(:,j),1)

            if Exx_mean(k)~=0
            E = [Exx_r(k,j) Exy_r(k,j); Exy_r(k,j) Eyy_r(k,j)];
            [vec,val] = eig(E);

            minVecs(:,k,j) = vec(:,1);
            maxVecs(:,k,j) = vec(:,2);
            minvals(k,j) = val(1,1);
            maxvals(k,j) = val(2,2);
            end
        end
        
end


minVecs_model = zeros(2,length(Exx_bs));
maxVecs_model = zeros(2,length(Exx_bs));
minvals_model = zeros(length(Exx_bs),1);
maxvals_model = zeros(length(Exx_bs),1);

for k=1:length(Exx_bs)

    
    E = [Exx_bs(k) Exy_bs(k); Exy_bs(k) Eyy_bs(k)];
    [vec,val] = eig(E);

    minVecs_model(:,k) = vec(:,1);
    maxVecs_model(:,k) = vec(:,2);
    minvals_model(k) = val(1,1);
    maxvals_model(k) = val(2,2);
    
end


mean_minvals = mean(minvals,2);
mean_maxvals = mean(maxvals,2);

vecscale = .10;
f=figure;
hold on

skip=5;

for k=1:length(Exx_bs)
    
    
        
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



        x0 = [0 vecscale*cos(t) 0]+obs_llh(k,1);
        z0 = [0 vecscale*sin(t) 0]+obs_llh(k,2);    
        if bigvals<0
             patch( x0, z0, 'r','FaceAlpha',.25 ,'EdgeColor','none');
        else
             patch( x0, z0, 'b','FaceAlpha',.25 ,'EdgeColor','none');
        end    

        x0 = [0 -vecscale*cos(t) 0]+obs_llh(k,1);
        z0 = [0 -vecscale*sin(t) 0]+obs_llh(k,2);

        if bigvals<0
           patch( x0, z0, 'r','FaceAlpha',.25 ,'EdgeColor','none');
        else
            patch( x0, z0, 'b','FaceAlpha',.25 ,'EdgeColor','none');
        end



        %plot smaller principal direction as a  line
        vx = [obs_llh(k,1)-mean_smallVecs(1)*vecscale*abs(smallvals/bigvals) obs_llh(k,1)+mean_smallVecs(1)*vecscale*abs(smallvals/bigvals)];
        vy = [obs_llh(k,2)-mean_smallVecs(2)*vecscale*abs(smallvals/bigvals) obs_llh(k,2)+mean_smallVecs(2)*vecscale*abs(smallvals/bigvals)];


        if smallvals<0
            plot(vx,vy,'r')
        else
            plot(vx,vy,'b')
        end

    end
   
    style(k) = (bigvals+smallvals)/(abs(bigvals)+abs(smallvals));
    
end



%plot model directions


for k=1:length(Exx_bs)
    
    
    %plot the larger magnitude principal direction
    
    if abs(minvals_model(k))>abs(maxvals_model(k))
        
        bigvals = minvals_model(k);
        smallvals = maxvals_model(k);
        
        bigVecs = minVecs_model(:,k);
        smallVecs = maxVecs_model(:,k);
        

                
    else
        
        bigvals = maxvals_model(k);
        smallvals = minvals_model(k);        
    
        bigVecs = maxVecs_model(:,k);
        smallVecs = minVecs_model(:,k);

    end
       
    
     model_angles(k) = atan(bigVecs(2,:)./bigVecs(1,:))*180/pi;
   
        if mod(k,skip)==0
   

            %plot larger principal direction 
            vx = [obs_llh(k,1)-bigVecs(1)*vecscale obs_llh(k,1)+bigVecs(1)*vecscale];
            vy = [obs_llh(k,2)-bigVecs(2)*vecscale obs_llh(k,2)+bigVecs(2)*vecscale];


            if bigvals<0
                plot(vx,vy,'r','linewidth',2)
            else
                plot(vx,vy,'b','linewidth',2)
            end


            %plot smaller principal direction 
            vx = [obs_llh(k,1)-smallVecs(1)*vecscale*abs(smallvals/bigvals) obs_llh(k,1)+smallVecs(1)*vecscale*abs(smallvals/bigvals)];
            vy = [obs_llh(k,2)-smallVecs(2)*vecscale*abs(smallvals/bigvals) obs_llh(k,2)+smallVecs(2)*vecscale*abs(smallvals/bigvals)];


            if smallvals<0
                plot(vx,vy,'r','linewidth',2)
            else
                plot(vx,vy,'b','linewidth',2)
            end

        end
        
         model_style(k) = (bigvals+smallvals)/(abs(bigvals)+abs(smallvals));

 end
    

plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('observed -- 95% wedges and thin lines,  model -- heavy lines') 
set(gca,'fontsize',15)

figure
subplot(121)
scatter(obs_llh(:,1),obs_llh(:,2),80,90-mean_angles,'fill')
colormap(pmkmp(128,'IsoAZ'))
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('observed azimith of larger principal strain rate')
set(gca,'fontsize',15)

subplot(122)
scatter(obs_llh(:,1),obs_llh(:,2),80,90-model_angles,'fill')
colormap(pmkmp(128,'IsoAZ'))
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('model azimith of larger principal strain rate')
set(gca,'fontsize',15)



figure
subplot(121)
scatter(obs_llh(:,1),obs_llh(:,2),80,-style,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('observed strain rate style (blue=normal, green=strike-slip, red=reverse')
set(gca,'fontsize',15)

subplot(122)
scatter(obs_llh(:,1),obs_llh(:,2),80,-model_style,'fill')
colormap(jet)
colorbar
hold on
plot(SegEnds_llh(:,[1 3])',SegEnds_llh(:,[2 4])','k','linewidth',1)
title('model strain rate style')
set(gca,'fontsize',15)


