%note: best to clear workspace because this scripts saves everything
clear all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SECTION

%origin of Cartensian coordinate system
%NOTE: this needs to be same origin used in generation of the mesh
%SegEnds cooredinates in FaultInfo
origin = [39 -121];




%% load fault info
% this needs to be a mat file that contains a cell array of FaultInfo
% structures for each fault in the model. For the k-th fault, FaultInfo{k}
% contains
%     FaultInfo{k}.SegEnds,  (Nx4 matrix of N fault segment coordinates (Cartesian) of
%       format [x1 x2 y1 y2] where (x1,y1) is first endpoint and (x2,y2)
%       is the other endpoint
%     FaultInfo{k}.FaultName, string, name of fault
%     FaultInfo{k}.FaultID, ID number associated with fault, scalar
%     FaultInfo{k}.Dip, dip angle (degrees), scalar
%     FaultInfo{k}.Dip_Bottom, dip angle (degrees) at bottom of fault, scalar
%     NOTE: Fault dip will vary linearlar from top (FaultInfo{k}.Dip) to
%           bottom (FaultInfo{k}.Dip_Bottom). FaultInfo{k}.Dip_Bottom defaults to 
%           uniform dip if not specified 
%     FaultInfo{k}.Rake, rake angle (degrees), scalar
%     FaultInfo{k}.PrefSlipRate, preffered slip rate (mm/yr), scalar  
%     FaultInfo{k}.UBSlipRate, upper bound on slip rate (mm/yr), scalar   
%     FaultInfo{k}.LBSlipRate, lower bound on slip rate (mm/yr), scalar    
%     FaultInfo{k}.STDSlipRate, standard deviation of slip rate (mm/yr), scalar   
%     NOTE: PrefSlipRate, UBSlipRate, LBSlipRate, or STDSlipRate can be set to
%          empty, [], if not known
%     FaultInfo{k}.RecurrenceTime, earthquake recurrence time, yrs, specify nan
%          if you don't want to compute earthquake cycles for the fault 
%          (i.e., use elastic halfspace solution), specify 'inf' to use the
%          default value specified below
%     FaultInfo{k}.TimeSinceEq, time since last last earthquake (yr),
%          specify 'inf' to use default setting specified beloe

%allow for variable rake? If true, two components of slip will be computed
%(twice the computation time
variable_rake = true;

%  average coseismic slip (meters) used to define default recurrence time 
coseismic_slip = 5;  %e.g., recurrence time = coseismic_slip*1000/PrefSlipRate

%default time since last earthqake as a fraction of the earthqake
%recurrence time
t_eq_frac = 0.2;

%name of fault info file
%faultinfo_name = './build_mesh/fault_info_NoCal_novisco.mat'; 
faultinfo_name = './build_mesh/fault_info_CSAF_novisco.mat'; 


%% mesh file name (mat file)
%meshfile =  './build_mesh/NoCal_mesh';
meshfile =  './build_mesh/CSAF_mesh';


%% slip patch discretization
nve_total = 3;  %number of elements down dip
pL=15; %nominal patch length (along strike)
DLock_max = 20;  %maximum depth of locking (depth to bottom of patches)
                 %Thtis can be a scalar (uniform for all faults) or a
                 %vector the length of number of fault segments for
                 %variable depths

%all variables will be save to floowing name of file (mat file) 
%build_filename = 'NoCal_bodyforce_moment_novisco';
build_filename = 'CSAF_bodyforce_moment_novisco';
                 

%load strainrate file (text file containing strain rate 'observations') 
%columns are: lon, lat, Exx, Exy, Eyy, Sig_Exx, Sig_Exy, Sig_Eyy
%If strain rates in units of 1/yr, set is_microstrain = false. If strain
%rates in micro-strain/yr, set is_microstrain = true
%
%NOTE: observations outside of model domain (meshed region) will be tossed
%out, or you can specify boundaries outside of which data is tossed out
strainrates_obs = './data/BforStrain_WUS_9_28_2022_smooth_trimesh.txt'; %note, micro-strain/yr
is_microstrain = true;

obs_limits = false;  %true or false, if true, data outside of limits (below) will be removed
obs_xmin = -300; obs_xmax = 0;
obs_ymin = -150; obs_ymax = 150;



%specify refinement for top row of slip patches (refined slip patch
%spacing, km)
refine = 1;  %shallow patch discretization (km)


% Strainrate Greens function are computed at centroids of mesh and interpolated to 
%   observation coordinates (in strainrates_obs)
% You may choose to keep strain rate GFs at centroids of mesh (might be large) for 
%   the purpose of plotting or you can discard to save space (true or
%   false)
save_strainrates_centroids = false;



%% viscoelastic cycle parameters
% number of terms in Hankel transform -- probably no need to change
Nterms = 100;  %number of terms in Hankel transfrom -- default is 100 and works well for many problems
             %where model domain is say <1000 km. May need to increase to 200-300  for
             %larger model domains. Aliasing of long wavelengths occurs if
             %too few terms. Use trial and error.

%layers
H1=DLock_max;  %thickness of elastic plate 
H2=50;  %depth to bottom of viscoelastic layer
tR1= 50; % relaxation time of viscoelastic layer (years)
tR2= 10;  %relaxation time of halfspace  (years)

%Poisson's ratio for distributed moment sources
nu = 0.25;


%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath backslip_tools
addpath tools
addpath visco_cycle


%% load everything

%load fault info
load(faultinfo_name)

%unwrap fault info and place in matrices/vectors
unwrap_faultinfo_visco



%load mesh
load(meshfile)

%load strainrates
A = load(strainrates_obs);

obs_llh = A(:,1:2);
xy_obs = llh2local(obs_llh',fliplr(origin))';

if is_microstrain
    scale=1e-6;
else
    scale=1;
end
Exx_mean = scale*A(:,3); 
Exy_mean = scale*A(:,4);
Eyy_mean = scale*A(:,5);

Exx_std = scale*A(:,6);
Exy_std = scale*A(:,7);
Eyy_std = scale*A(:,8);

%% toss out strain rate obs outside of model domain
minx = min(tri_centroids(:,1));
maxx = max(tri_centroids(:,1));
miny = min(tri_centroids(:,2));
maxy = max(tri_centroids(:,2));


outside = xy_obs(:,1)<minx |  xy_obs(:,1)>maxx |  xy_obs(:,2)<miny | xy_obs(:,2)>maxy;
xy_obs(outside,:) = [];
obs_llh(outside,:) = [];
Exx_mean(outside) = [];
Exy_mean(outside) = [];
Eyy_mean(outside) = [];
Exx_std(outside) = [];
Exy_std(outside) = [];
Eyy_std(outside) = [];


%specified boundaries for data
if obs_limits
    outside = xy_obs(:,1)<obs_xmin |  xy_obs(:,1)>obs_xmax |  xy_obs(:,2)<obs_ymin | xy_obs(:,2)>obs_ymax;
    xy_obs(outside,:) = [];
    obs_llh(outside,:) = [];
    Exx_mean(outside) = [];
    Exy_mean(outside) = [];
    Eyy_mean(outside) = [];
    Exx_std(outside) = [];
    Exy_std(outside) = [];
    Eyy_std(outside) = [];
end


%% extract triangular mesh info 
patch_stuff=make_triangular_patch_stuff(tri,[nodes 0*nodes(:,1)]);
tri_areas = patch_stuff.area_faces;
tri_centroids = patch_stuff.centroids_faces;

    
%% build slip patches
Dips = Seg_Dips;
if length(DLock_max)==1
    DLock_max = DLock_max*ones(size(Dips));
end

make_patches_backslip_toprow_variableH
make_patches_backslip_deep_variableH




%%Build Greens functions (piece-wise linear slip gradients)

%top row patches (piece-wise linear slip gradients)
[G1Exx_top,G2Exx_top,G1Exy_top,G2Exy_top,...
    G1Eyy_top,G2Eyy_top] = make_strainG_piecewise(pm_top,rake_top_patches,tri_centroids(:,1:2),refine);
%combine matrices above to build G matrices
[GExx_top_fine,GExy_top_fine,GEyy_top_fine]=...
    build_G_backslip_function(G1Exx_top,G2Exx_top,G1Exy_top,G2Exy_top,G1Eyy_top,G2Eyy_top,pm_top_seg_num);
%clear variables no longer needed
clear G1Exx_top G2Exx_top G1Exy_top G2Exy_top G1Eyy_top G2Eyy_top


if variable_rake
    
    %add 90 to rake to get slip component perpedicular to specificied rake
    %diretion
    
    %top row patches (piece-wise linear slip gradients)
    [G1Exx_top,G2Exx_top,G1Exy_top,G2Exy_top,...
        G1Eyy_top,G2Eyy_top] = make_strainG_piecewise(pm_top,rake_top_patches+90,tri_centroids(:,1:2),refine);
    %combine matrices above to build G matrices
    [GExx_top_fine_perp,GExy_top_fine_perp,GEyy_top_fine_perp]=...
        build_G_backslip_function(G1Exx_top,G2Exx_top,G1Exy_top,G2Exy_top,G1Eyy_top,G2Eyy_top,pm_top_seg_num);
    %clear variables no longer needed
    clear G1Exx_top G2Exx_top G1Exy_top G2Exy_top G1Eyy_top G2Eyy_top
    
end

%G matrices for deeper patches (below top row)
[GExx_fine,GExy_fine,GEyy_fine] = make_strainG(pm,rake_patches,tri_centroids(:,1:2));

if variable_rake
    
    %add 90 to rake to get slip component perpedicular to specificied rake
    %diretion
   
    [GExx_fine_perp,GExy_fine_perp,GEyy_fine_perp] = make_strainG(pm,rake_patches+90,tri_centroids(:,1:2));

end

%interpolate to observation coordinates
GExx_top = zeros(size(xy_obs,1),size(pm_top,1));
GExy_top = zeros(size(xy_obs,1),size(pm_top,1));
GEyy_top = zeros(size(xy_obs,1),size(pm_top,1));
 
GExx = zeros(size(xy_obs,1),size(pm,1));
GExy = zeros(size(xy_obs,1),size(pm,1));
GEyy = zeros(size(xy_obs,1),size(pm,1));

if variable_rake

    GExx_top_perp = zeros(size(xy_obs,1),size(pm_top,1));
    GExy_top_perp = zeros(size(xy_obs,1),size(pm_top,1));
    GEyy_top_perp = zeros(size(xy_obs,1),size(pm_top,1));

    GExx_perp = zeros(size(xy_obs,1),size(pm,1));
    GExy_perp = zeros(size(xy_obs,1),size(pm,1));
    GEyy_perp = zeros(size(xy_obs,1),size(pm,1));

end

for k=1:size(GExx_top,2)

    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_top_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_top_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GExx_top(:,k) = grid;

    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_top_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_top_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GExy_top(:,k) = grid;
    
    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_top_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_top_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GEyy_top(:,k) = grid;
   
    if variable_rake
        
        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_top_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_top_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GExx_top_perp(:,k) = grid;

        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_top_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_top_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GExy_top_perp(:,k) = grid;

        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_top_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_top_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GEyy_top_perp(:,k) = grid;
   
    end

    disp(['Completed ' num2str(k/size(GExx_top,2)*100) '% of top row interpolations'])
end

for k=1:size(GExx,2)
   
    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GExx(:,k) = grid;

    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GExy(:,k) = grid;

    grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_fine(:,k),xy_obs(:,1),xy_obs(:,2));
    grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_fine(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
    GEyy(:,k) = grid;

     if variable_rake
         
        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExx_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GExx_perp(:,k) = grid;

        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GExy_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GExy_perp(:,k) = grid;

        grid = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_fine_perp(:,k),xy_obs(:,1),xy_obs(:,2));
        grid(isnan(grid)) = griddata(tri_centroids(:,1),tri_centroids(:,2),GEyy_fine_perp(:,k),xy_obs(isnan(grid),1),xy_obs(isnan(grid),2),'nearest');
        GEyy_perp(:,k) = grid;
    
     end
     
    disp(['Completed ' num2str(k/size(GExx,2)*100) '% of deeper row interpolations'])

end




%viscoelastic GFs
%Top row

Ts_top = zeros(size(pm_top_seg_id));
teqs_top = nan(size(pm_top_seg_id));
pref_rate_pm_top = nan(size(pm_top_seg_id));

for k=1:length(Seg_FaultID)
    
    ind = pm_top_seg_id==Seg_FaultID(k);
    Ts_top(ind) = Seg_RecurrenceTime(k);
    teqs_top(ind) = Seg_TimeSinceEq(k);
    pref_rate_pm_top(ind) = Seg_GeoPrefRate(k);
    
end

ind = isinf(Ts_top);
Ts_top(ind) = coseismic_slip*1000./pref_rate_pm_top(ind);
ind = isinf(teqs_top);
teqs_top(ind) = t_eq_frac*Ts_top(ind);

%viscoelastic
[GExx_top_cycle,GExy_top_cycle,GEyy_top_cycle]=make_StrainG_cycles(pm_top,rake_top_patches,...
    xy_obs,H1,H2,tR1,tR2,Ts_top,teqs_top,Nterms);

if variable_rake
    
    %add 90 to rake to get slip component perpedicular to specificied rake
    %diretion
    [GExx_top_cycle_perp,GExy_top_cycle_perp,GEyy_top_cycle_perp]=make_StrainG_cycles(pm_top,rake_top_patches+90,...
        xy_obs,H1,H2,tR1,tR2,Ts_top,teqs_top,Nterms);
end

%bottom rows
Ts = zeros(size(pm_seg_id));
teqs = nan(size(pm_seg_id));
pref_rate_pm = nan(size(pm_seg_id));

for k=1:length(Seg_FaultID)
    
    ind = pm_seg_id==Seg_FaultID(k);
    Ts(ind) = Seg_RecurrenceTime(k);
    teqs(ind) = Seg_TimeSinceEq(k);
    pref_rate_pm(ind) = Seg_GeoPrefRate(k);
    
end

ind = isinf(Ts);
Ts(ind) = coseismic_slip*1000./pref_rate_pm(ind);
ind = isinf(teqs);
teqs(ind) = t_eq_frac*Ts(ind);

%viscoelastic
[GExx_cycle,GExy_cycle,GEyy_cycle]=make_StrainG_cycles(pm,rake_patches,...
    xy_obs,H1,H2,tR1,tR2,Ts,teqs,Nterms);

if variable_rake
    
    %add 90 to rake to get slip component perpedicular to specificied rake
    %diretion
    [GExx_cycle_perp,GExy_cycle_perp,GEyy_cycle_perp]=make_StrainG_cycles(pm,rake_patches+90,...
    xy_obs,H1,H2,tR1,tR2,Ts,teqs,Nterms);

end

GExx_top_elastic = -GExx_top; %negative converts to backslip
GExy_top_elastic = -GExy_top;
GEyy_top_elastic = -GEyy_top;
GExx_elastic = -GExx;
GExy_elastic = -GExy;
GEyy_elastic = -GEyy;

if variable_rake
    GExx_top_elastic_perp = -GExx_top_perp; %negative converts to backslip
    GExy_top_elastic_perp = -GExy_top_perp;
    GEyy_top_elastic_perp = -GEyy_top_perp;
    GExx_elastic_perp = -GExx_perp;
    GExy_elastic_perp = -GExy_perp;
    GEyy_elastic_perp = -GEyy_perp;
end


%build 2d moment source GFs (distributed moment sources)
[Gexx_mom,Gexy_mom,Geyy_mom] = buildG_MomentSource_2d(xy_obs,tri, nodes,tri_centroids,nu);

%only keep the strain rates on the fine triangular based mesh if requested
%these are large matrices 
if ~save_strainrates_centroids
   clear GExx_top_fine GExy_top_fine GEyy_top_fine GExx_fine GExy_fine GEyy_fine
   if variable_rake
        clear GExx_top_fine_perp GExy_top_fine_perp GEyy_top_fine_perp GExx_fine_perp GExy_fine_perp GEyy_fine_perp
   end
   
end

clear variable_rake
%save everything to a mat file
save(build_filename)