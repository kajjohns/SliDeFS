
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUT SECTION

%name of file generated from build_backslip_GreenFunctions
build_filename = 'CSAF';

%regular grid? True or false. If true, provide grid spacing in degrees
%below. If grid is not regular, specify false
regular_grid = false;
dl =0.2;  %grid spacing (degrees, needed if above is true)

%% choose type of inversion

%allow for variable rake? If true, two components of slip will be computed
%(twice the computation time, if false, rake is fixed 
variable_rake = false;

%weight placed on minimizing deviation from rake (if true above)
rake_weight = 100;

%use elastic-only solution (ignore viscoelastic cycle GFs)
use_elastic = true;

%Gaussian slip deficit rate prior?
use_Gaussian_prior = false;  %true/false, use Guassian slip rate prior
                            %NOTE: if true, pref_slip_rate must be defined 
                            %if true, std_slip_rate must be defined or set use_bounds_for_std = true;
use_bounds_for_std = true;  %true/false (only relevant if above is true), if true std of slip rate is (ub-lb)*std_scale
std_scale = .5;

%use bounds on slip deficit rate?   
use_upper_bounds = true;  %true to use upper bounds
                          %NOTE: if ture, ub_slip_rate must be defined
use_lower_bounds = false;  %true to use lower bounds
                           %if false, lower bounds will be zero
                           %NOTE: if true, lb_slip_rate must be defined
scale_upper_bounds = 10;  %optional scaling of upper bounds, =1 for no scaling
scale_lower_bounds = 1;   %optional scaling of lower bounds, =1 for no scaling

%smoothing weight (weight placed on keeping slip deficit rate smooth vs.
%fitting data
weight_smooth = .001;  
                         

%%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale weights to order of magntiude of strain rates
rake_weight = 10^6*rake_weight;
weight_smooth = 10^6*weight_smooth;

addpath backslip_tools
addpath tools

%load Greens Functions
load(build_filename)

%make smoothing operator
make_pm_smoothing

%check for variable rake GFs
if variable_rake
    if ~exist('GExx_top_cycle_perp')
        disp('Variable rake GFs were not computed. Script terminated.')
        return
    end
end


if use_elastic
   
    GExx_top = GExx_top_elastic;
    GExy_top = GExy_top_elastic;
    GEyy_top = GEyy_top_elastic;

    GExx = GExx_elastic;
    GExy = GExy_elastic;
    GEyy = GEyy_elastic;

    if variable_rake
        GExx_top_perp = GExx_top_elastic_perp;
        GExy_top_perp = GExy_top_elastic_perp;
        GEyy_top_perp = GEyy_top_elastic_perp;

        GExx_perp = GExx_elastic_perp;
        GExy_perp = GExy_elastic_perp;
        GEyy_perp = GEyy_elastic_perp;
    end
    
else
    
    GExx_top = GExx_top_elastic + 1000*GExx_top_cycle*10^-6; %convert from micro-strain/yr, convert from m/yr to mm/yr
    GExy_top = GExy_top_elastic + 1000*GExy_top_cycle*10^-6;
    GEyy_top = GEyy_top_elastic + 1000*GEyy_top_cycle*10^-6;

    GExx = GExx_elastic + 1000*GExx_cycle*10^-6;
    GExy = GExy_elastic + 1000*GExy_cycle*10^-6;
    GEyy = GEyy_elastic + 1000*GEyy_cycle*10^-6;

    if variable_rake
        
        GExx_top_perp = GExx_top_elastic + 1000*GExx_top_cycle_perp*10^-6; %convert from micro-strain/yr, convert from m/yr to mm/yr
        GExy_top_perp = GExy_top_elastic + 1000*GExy_top_cycle_perp*10^-6;
        GEyy_top_perp = GEyy_top_elastic + 1000*GEyy_top_cycle_perp*10^-6;

        GExx_perp = GExx_elastic_perp + 1000*GExx_cycle_perp*10^-6;
        GExy_perp = GExy_elastic_perp + 1000*GExy_cycle_perp*10^-6;
        GEyy_perp = GEyy_elastic_perp + 1000*GEyy_cycle_perp*10^-6;
    end
    

end



G = [GExx_top GExx;...
      GExy_top GExy; ...
      GEyy_top GEyy];

if variable_rake
    
    G_perp = [GExx_top_perp GExx_perp;...
    GExy_top_perp GExy_perp; ...
    GEyy_top_perp GEyy_perp];
  
    G = [G G_perp];
    
end

    
%bounds
if use_upper_bounds
    UB = [ub_top_patches; ub_patches]/1000;  %convert to m/yr
    if variable_rake; UB = [UB; UB]; end
else
    UB = [];
end

if use_lower_bounds
    LB = [lb_top_patches; lb_patches]/1000;  %convert to m/yr
    if variable_rake; LB = [LB; LB]; end
else
    LB = [0*ub_top_patches; 0*ub_patches];
    if variable_rake; LB = [LB; LB]; end
end

UB = UB*scale_upper_bounds;
LB = LB*scale_lower_bounds;

if use_Gaussian_prior
    mean_s_prior = [pref_top_patches; pref_patches]/1000;  %convert to m/yr;
    if use_bounds_for_std 
        std_s_prior = std_scale*([ub_top_patches; ub_patches]-[lb_top_patches; lb_patches])/1000;  %convert to m/yr;
        if min(std_s_prior)<=0
            std_s_prior(std_s_prior<=0)=10/1000;  %convert to m/yr;
            disp('Warning: zero or negative standard deviation set to 10 mm/yr')
        end

    else
        std_s_prior = [std_top_patches; std_patches]/1000;  %convert to m/yr
    end
end


%conduct weighted least squares    
Esig = [Exx_std;Exy_std;Eyy_std];
d = [Exx_mean;Exy_mean;Eyy_mean];

if use_Gaussian_prior
    Gprior = diag(1./std_s_prior);
    dd = [d./Esig; mean_s_prior./std_s_prior;zeros(size(pm_smooth,1),1)];
    if variable_rake
        GG = [G./repmat(Esig,1,size(G,2));Gprior;weight_smooth*[pm_smooth pm_smooth]];
    else
        GG = [G./repmat(Esig,1,size(G,2));Gprior;weight_smooth*pm_smooth];
    end
else
    dd = [d./Esig;zeros(size(pm_smooth,1),1)];
     if variable_rake
        GG = [G./repmat(Esig,1,size(G,2));weight_smooth*[pm_smooth pm_smooth]];
     else
         GG = [G./repmat(Esig,1,size(G,2));weight_smooth*pm_smooth];
     end   
end

if variable_rake
   %damp rake-perp component of slip 
   N = size(pm,1)+size(pm_top,1);
   Gperp = [zeros(N,N) rake_weight*eye(N,N)];
   GG = [GG; Gperp];
   dd = [dd; zeros(N,1)];
end


%bounded least squares result
mhat = lsqlin(GG,dd,[],[],[],[],LB,UB);


%predicted strain rates
dhat = G*mhat;

%reduced chi-squared
chi2 = (d./Esig-dhat./Esig)'*(d./Esig-dhat./Esig)/length(d)

data_Var = (d)'*(d);
resid_Var = (d-dhat)'*(d-dhat);

var_reduction = 1 - resid_Var/data_Var

    
    
    
    


% compute moments from backslip
A_top = pm_top(:,1).*pm_top(:,2)*10^6;
A_bot = pm(:,1).*pm(:,2)*10^6;

  
        
Exx_bs = dhat(1:end/3);
Exy_bs = dhat(1+end/3:2*end/3);
Eyy_bs = dhat(1+2*end/3:end);
  
    
Exx_offfault = Exx_mean - Exx_bs;
Exy_offfault = Exy_mean - Exy_bs;
Eyy_offfault = Eyy_mean - Eyy_bs;

max_shear = 10^6*sqrt((Exx_mean-Eyy_mean).^2 + Exy_mean.^2); %convert to micro-strain/yr
max_shear_bs = 10^6*sqrt((Exx_bs-Eyy_bs).^2 + Exy_bs.^2);

dilatation = 10^6*(Exx_mean + Eyy_mean);
dilatation_bs = 10^6*(Exx_bs + Eyy_bs);

max_shear_offfault = 10^6*sqrt(((Exx_mean-Exx_bs)-(Eyy_mean-Eyy_bs)).^2 + (Exy_mean-Exy_bs).^2);
dilatation_offfault = 10^6*((Exx_mean-Exx_bs)+(Eyy_mean-Eyy_bs));

%compute moments
H=15;
mu=30e9;
% 
% 
% Kostrov_Mo_total =    Kostrov(xy_obs(:,1),xy_obs(:,2),Exx_mean,Exy_mean,Eyy_mean,H,mu);
% Kostrov_Mo_bs =       Kostrov(xy_obs(:,1),xy_obs(:,2),Exx_bs,Exy_bs,Eyy_bs,H,mu);
% Kostrov_Mo_offfault = Kostrov(xy_obs(:,1),xy_obs(:,2),Exx_offfault,Exy_offfault,Eyy_offfault,H,mu);
% 
%     

if variable_rake
    mhat_rake = mhat(1:end/2); %slip in rake direction
    mhat_perp = mhat(1+end/2:end); %slip in rake-perpendicular direction
else
    mhat_rake = mhat;
end

bs_top_rake = mhat_rake(1:size(pm_top,1));
bs_bot_rake = mhat_rake(1+size(pm_top,1):end);

if variable_rake
    bs_top_perp = mhat_perp(1:size(pm_top,1));
    bs_bot_perp = mhat_perp(1+size(pm_top,1):end);
end

if variable_rake
    bs_top = sqrt(bs_top_rake.^2 + bs_top_perp.^2);
    bs_bot = sqrt(bs_bot_rake.^2 + bs_bot_perp.^2);
else
    bs_top = bs_top_rake;
    bs_bot = bs_bot_rake;
end


if variable_rake
    %slip vector in plane -- compute rake
    v_r_top_x = bs_top_rake.*cos(rake_top_patches*pi/180);
    v_r_top_y = bs_top_rake.*sin(rake_top_patches*pi/180);
    v_p_top_x = bs_top_perp.*cos((rake_top_patches+90)*pi/180);
    v_p_top_y = bs_top_perp.*sin((rake_top_patches+90)*pi/180);
    v_tot_top_x = v_r_top_x + v_p_top_x;
    v_tot_top_y = v_r_top_y + v_p_top_y;
    rake_top = atan2(v_tot_top_y,v_tot_top_x)*180/pi;

    v_r_bot_x = bs_bot_rake.*cos(rake_patches*pi/180);
    v_r_bot_y = bs_bot_rake.*sin(rake_patches*pi/180);
    v_p_bot_x = bs_bot_perp.*cos((rake_patches+90)*pi/180);
    v_p_bot_y = bs_bot_perp.*sin((rake_patches+90)*pi/180);
    v_tot_bot_x = v_r_bot_x + v_p_bot_x;
    v_tot_bot_y = v_r_bot_y + v_p_bot_y;
    rake_bot = atan2(v_tot_bot_y,v_tot_bot_x)*180/pi;

else
    
    rake_top = rake_top_patches;
    rake_bot = rake_patches;

end


Mo_fault_rec = mu*sum(A_top.*bs_top)+mu*sum(A_bot.*bs_bot);
Mo_fault_rec_bot = mu*sum(A_bot.*bs_bot);


SegEnds_llh1 = local2llh(SegEnds(:,1:2)',fliplr(origin))';
SegEnds_llh2 = local2llh(SegEnds(:,3:4)',fliplr(origin))';
SegEnds_llh = [SegEnds_llh1 SegEnds_llh2];
