
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUT SECTION

%name of file generated from build_backslip_GreenFunctions
build_filename = 'NoCal_bodyforce_moment_novisco';

%regular grid? True or false. If true, provide grid spacing in degrees
%below. If grid is not regular, specify false
regular_grid = false;
dl =0.2;  %grid spacing (degrees, needed if above is true)

%% choose type of inversion

%allow for variable rake? If true, two components of slip will be computed
%(twice the computation time, if false, rake is fixed 
variable_rake = true;

%weight placed on minimizing deviation from rake (if true above)
weight_rake = 1;

%use elastic-only solution (ignore viscoelastic cycle GFs)
use_elastic = false;

%Gaussian slip deficit rate prior?
use_Gaussian_prior = true;  %true/false, use Guassian slip rate prior
                            %NOTE: if true, pref_slip_rate must be defined 
                            %if true, std_slip_rate must be defined or set use_bounds_for_std = true;
use_bounds_for_std = true;  %true/false (only relevant if above is true), if true std of slip rate is (ub-lb)*std_scale
std_scale = .1;

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
weight_smooth = 1;  

%include distributed moment sources?
include_moment = true;

%damping weight for moment forces (weight to keep moment forces small)
weight_moment = 1;  


%%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale weights to improve relative weighting
weight_moment = 10^3*weight_moment;
weight_smooth = 10^3*weight_smooth;
weight_rake = 10^3*weight_rake;


addpath backslip_tools
addpath tools

%load Greens Functions
load(build_filename)

%make smoothing operator
make_pm_smoothing


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

if include_moment

    %moment sources
    Gmom = [Gexx_mom.m11 Gexx_mom.m12 Gexx_mom.m22];
    Gmom = [Gmom; [Gexy_mom.m11 Gexy_mom.m12 Gexy_mom.m22]];
    Gmom = [Gmom; [Geyy_mom.m11 Geyy_mom.m12 Geyy_mom.m22]];
    Gmom = Gmom*10^-3;  %improve scaling
    
    G = [G Gmom];  

else
    Gmom = [];
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
    if variable_rake; LB = [LB; -LB]; end
else
    LB = [0*ub_top_patches; 0*ub_patches];
    if variable_rake; LB = [LB; -UB]; end
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


%backslip weighting matrix, W
if variable_rake
        W = [weight_smooth*[pm_smooth pm_smooth] zeros(size(pm_smooth,1),size(Gmom,2))];
else
        W = [weight_smooth*pm_smooth zeros(size(pm_smooth,1),size(Gmom,2))];
end

if use_Gaussian_prior

    if variable_rake
    
        Gprior = [diag(1./std_s_prior) zeros(length(std_s_prior))];
        Gprior = [Gprior zeros(size(Gprior,1),size(Gmom,2))];  %zeros for moment source terms
        dd = [d./Esig; mean_s_prior./std_s_prior;zeros(size(pm_smooth,1),1)];
        GG = [G./repmat(Esig,1,size(G,2));Gprior;W];

    else

        Gprior = diag(1./std_s_prior);
        Gprior = [Gprior zeros(size(Gprior,1),size(Gmom,2))];  %zeros for moment source terms
        dd = [d./Esig; mean_s_prior./std_s_prior;zeros(size(pm_smooth,1),1)];
        GG = [G./repmat(Esig,1,size(G,2));Gprior;W];

    end

else

    dd = [d./Esig;zeros(size(pm_smooth,1),1)];
    GG = [G./repmat(Esig,1,size(G,2));W];

end

%damp rake-perp component of slip
N = size(pm,1)+size(pm_top,1); %number of slip patches
if variable_rake
   Gdamp_perp = [zeros(N,N) weight_rake*eye(N,N) zeros(N,size(Gmom,2))];
   GG = [GG; Gdamp_perp];
   dd = [dd; zeros(N,1)];
end

%damp moment forces
M = size(Gmom,2);
if variable_rake
    Gdamp = [zeros(M,2*N) weight_moment*eye(M)];
else
    Gdamp = [zeros(M,N) weight_moment*eye(M)];
end

GG = [GG; Gdamp];
dd = [dd; zeros(size(Gdamp,1),1)];


%bounded least squares result
warning off  %turn off warnings about bounds being set to +/- inf which is totally fine, if not a bit sloppy
options = optimoptions(@lsqlin,'Display','iter');
display('lsqlin iteration status:')
mhat = lsqlin(GG,dd,[],[],[],[],LB,UB,[],options);
warning on

mhat_bs = mhat(1:end-M);
mhat_mom = mhat(end-M+1:end);

%predicted strain rates
dhat = G*mhat;

%reduced chi-squared
chi2 = (d./Esig-dhat./Esig)'*(d./Esig-dhat./Esig)/length(d)

%variance reduction
data_Var = (d)'*(d);
resid_Var = (d-dhat)'*(d-dhat);
var_reduction = 1 - resid_Var/data_Var

   

% compute moments from backslip
A_top = pm_top(:,1).*pm_top(:,2)*10^6;
A_bot = pm(:,1).*pm(:,2)*10^6;

Exx_total = dhat(1:end/3);
Exy_total = dhat(1+end/3:2*end/3);
Eyy_total = dhat(1+2*end/3:end);

if variable_rake
    dhat_bs = G(:,1:2*N)*mhat_bs;
else
    dhat_bs = G(:,1:N)*mhat_bs;
end

%strain rates from backslip
Exx_bs = dhat_bs(1:end/3);
Exy_bs = dhat_bs(1+end/3:2*end/3);
Eyy_bs = dhat_bs(1+2*end/3:end);
 
%strain rates from distributed sources
dhat_mom = Gmom*mhat_mom;  
Exx_bf = dhat_mom(1:end/3);
Exy_bf = dhat_mom(1+end/3:2*end/3);
Eyy_bf = dhat_mom(1+2*end/3:end);


max_shear = 10^6*sqrt((Exx_mean-Eyy_mean).^2 + Exy_mean.^2); %observed, convert to micro-strain/yr
max_shear_bs = 10^6*sqrt((Exx_bs-Eyy_bs).^2 + Exy_bs.^2);
max_shear_bf = 10^6*sqrt((Exx_bf-Eyy_bf).^2 + Exy_bf.^2);
max_shear_total = 10^6*sqrt((Exx_total-Eyy_total).^2 + Exy_total.^2);

dilatation = 10^6*(Exx_mean + Eyy_mean);
dilatation_bs = 10^6*(Exx_bs + Eyy_bs);
dilatation_bf = 10^6*(Exx_bf + Eyy_bf);
dilatation_total = 10^6*(Exx_total + Eyy_total);

max_shear_residual = 10^6*sqrt(((Exx_mean-Exx_total)-(Eyy_mean-Eyy_total)).^2 + (Exy_mean-Exy_total).^2);
dilatation_residual = 10^6*((Exx_mean-Exx_total)+(Eyy_mean-Eyy_total));



if variable_rake
    mhat_rake = mhat_bs(1:end/2); %slip in rake direction
    mhat_perp = mhat_bs(1+end/2:end); %slip in rake-perpendicular direction
else
    mhat_rake = mhat_bs;
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

%% compute on-fault and off-fault moment
H=20; %effective elastic thickness
mu=30e9; %shear modulus
 

%on fault
%note about units; bs i m/yr, Areas in meters, so Mo is N*m/yr
Mo_fault = mu*sum(A_top.*bs_top)+mu*sum(A_bot.*bs_bot) + mu*sum(A_bot.*bs_bot);

%off fault (distributed moments)

%Here we want to convert the weights (mhat_mom) on the moment sources to units of
%Moment [Force*distance, N*m]. The weights are moment normalized by plate
%thickness, H, and elastic shear modulus, mu
%
%two considerations in converting sources to moments in units of N*m
%1. force-couple GFs (Gmom matrix) are in units of 1/L^2 where Length is in km (need to
%convert to meters). But note that Gmom was already rescaled (for improve scaling
%for lsq inversion) as Gmom = Gmom*10^-3. So another factor of 10^-3 needs
%to be applied to convert GFs to units of 1/m^2
%2. After considering conversion of Gmom in step 1 (i.e., multiply weights by 1e3), weights on 
%moment terms are in units of moment normalized by mu*H.  


%convert mhat_mom to moments
mu = 3.0000e+10;
H = 25000; %average crustal thickness, meters
conv = 1e3*mu*H;  %conversion factor (see notes above)
m11 = mhat_mom(1:end/3)*conv; %dipole, force couple
m12_m21 = mhat_mom(1+end/3:2*end/3)*conv; %this is actually double couple m12+m21
m22 = mhat_mom(1+2*end/3:end)*conv; %dipole, force couple

%compute deviatoric moment tensor to get moment (Mo) of double-couple representions
for j=1:length(m11)
    E = [m11(j)-(m11(j)+m22(j))/2 .5*m12_m21(j); .5*m12_m21(j) m22(j)-(m11(j)+m22(j))/2] ;
    [vec,val] = eig(E);
    Mo(j) = abs(val(1,1));
end

Mo_off = sum(Mo);  %off-fault moment computed as the total double couple moment


%print moments
disp(['On-fault moment accumulation rate: ' num2str(Mo_fault) ' Nm/yr'])
disp(['Off-fault moment rate (double couples): ' num2str(Mo_off) ' Nm/yr'])
disp(['Total (on+ff) moment rate : ' num2str(Mo_fault + Mo_off) ' Nm/yr'])
disp(['Percent of total due to on-fault: ' num2str(Mo_fault/(Mo_fault + Mo_off)*100) '%'])
disp(['Percent of total due to off-fault: ' num2str(Mo_off/(Mo_fault + Mo_off)*100) '%'])

SegEnds_llh1 = local2llh(SegEnds(:,1:2)',fliplr(origin))';
SegEnds_llh2 = local2llh(SegEnds(:,3:4)',fliplr(origin))';
SegEnds_llh = [SegEnds_llh1 SegEnds_llh2];

nodes_llh = local2llh(nodes',fliplr(origin))';
