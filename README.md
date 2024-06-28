[![DOI](https://zenodo.org/badge/616735576.svg)](https://zenodo.org/badge/latestdoi/616735576)

# SliDeFS
SliDeFS (Slip Deficit From Strain) is a method for inverting surface strain rate components for slip deficit rate on faults as described in Johnson et al. (2024, see citation at bottom). The most recent release of SliDeFS inverts for slip deficit rates on faults and distributed moment sources. Faults are discretized into rectangular slip patches and the Okada (1985) solution is used to relate slip deficit rate (backslip) to surface strain rates (Exx, Exy, Eyy). Off-fault moment sources are represented with analytical solutions for force couples in a thin elastic plate (plane stress elasticity). The inversion is a bounded least squares inversion with optional spatial smoothing of slip deficit rate and damping of moment sources. Gaussian or truncated Gaussian priors on slip deficit rate can be included. Strain rates are computed on a triangular mesh constructed from surface traces of the fault model to minimize strain singularity artifacts. Details are provided in Johnson et al. (2022). Moment sources are computed at nodes of the same triangular mesh used for the slip deficit rate calculations. This mesh can be quite dense in areas of fault complexity. It would not be difficult to modify the code to utilize a different, simpler, mesh for the moment sources, but this has not been coded into the current release. 

A brief listing of the workflow steps is provided here followed by more detailed explanations below. 

1. Prepare fault geometry and create FaultInfo mat file. Prepare strain rate observation file (text file) with columns
lon, lat, Exx, Exy, Eyy, Sig_Exx, Sig_Exy, Sig_Eyy
where strain rates and standard deviations are in units of 1/yr or   micro-strain/yr (units are specified with a flag in the code)
2. In the build_mesh folder, run build_mesh.m to create the triangular mesh for strain rate Greens functions. 

The remaining steps are all conducted in the main SliDeFS folder. 

3. Run build_backslip_GreensFunctions.m.  Will take a few minutes for relatively small problem (100-200 segments) up to several hours for a large problem (thousands of segments).
4. Run invert_strainrate_for_backslip.m. Will take a few seconds for small problems and up to ~15-20 minutes for large problems. 
5. Run plot_best_fitting.m to view inversion results. 
6. (optional) Run plot_principal_directions.m to view observed and model principal strain rate directions
7. (optional) Run sample_posterior.m to compute posterior distribution of slip deficit rate.  (minutes to hours depending on size of problem)
8. (optional – requires step 7) Run plot_slip_rates_along_faults.m to visualize depth-averaged slip deficit rates along specified faults. 
9. (optional – requires step 7) Run plot_moment_sources.m to visualize the spatial distribution of principal moment rate sources and uncertainties. 







More detailed description of workflow. 

1. Prepare the fault geometry and create FaultInfo file. I have included an example, fault_info_CSAF.mat, for central San Andreas Fault. The method assumes rectangular, planar fault sections defined with a strike and dip using the right-hand rule -- thumb with palm down points in direction of strike (thumb points from first endpoint of segment towards second endpoint)  and index finger in direction of positive dip angle (between 0 and 90 degrees). You need to build the FaultInfo mat file that contains a cell array of fault information as shown below. Each entry of the array is for a fault comprised of smaller segments (often many). Rake is defined as the angle in the fault plane between the slip vector and the strike vector. Zero degree rake is left-lateral, +/-180 is right lateral, 90 is reverse and -90 is normal sense of slip.      

% this needs to be a mat file that contains a cell array of FaultInfo
% structures for each fault in the model. For the k-th fault, FaultInfo{k}
% contains
%     FaultInfo{k}.SegEnds,  (Nx4 matrix of N fault segment coordinates (Cartesian) of
%       format [x1 x2 y1 y2] where (x1,y1) is first endpoint and (x2,y2)
%       is the other endpoint
%     FaultInfo{k}.FaultName, string, name of fault
%     FaultInfo{k}.FaultID, ID number associated with fault, scalar
%     FaultInfo{k}.Dip, dip angle (degrees, positive), scalar
%     FaultInfo{k}.Rake, rake angle (degrees), scalar
%     FaultInfo{k}.PrefSlipRate, preferred slip rate (mm/yr), scalar  
%     FaultInfo{k}.UBSlipRate, upper bound on slip rate (mm/yr), scalar   
%     FaultInfo{k}.LBSlipRate, lower bound on slip rate (mm/yr), scalar    
%     FaultInfo{k}.STDSlipRate, standard deviation of slip rate (mm/yr), scalar   
%     NOTE: PrefSlipRate, UBSlipRate, LBSlipRate, or STDSlipRate can be set to
%          empty, [], if not known
%     FaultInfo{k}.RecurrenceTime, earthquake recurrence time, yrs, specify nan
%          if you don't want to compute earthquake cycles for the fault 
%          (i.e., use elastic halfspace solution), specify 'inf' to use the
%          default value 
%     FaultInfo{k}.TimeSinceEq, time since last last earthquake (yr),
%          specify 'inf' to use default setting 


2. Create the triangular mesh for strain rate Greens functions using build_mesh.m. Various parameters are specified in INPUT SECTION. This script requires the FaultInfo file because the mesh depends on fault segment geometry. Specify the nominal slip patch length (along strike) – this needs to be the same length use to make the slip patches in build_backslip_GreensFunctions.m. hmax is the largest allowable patch size in the mesh (side length of triangle). The mesh geometry parameters are saved in mat file named by string variable ‘savename’.

3. Build Greens Functions. Various parameters and options are set in the INPUT SECTION: 
(a) The origin (lat, lon) is used to convert to Cartesian coordinate system and must be same origin used to for the Cartesian geometry specified in FaultInfo. 
(b) variable_rake is ‘true’ if you want to allow estimated slip deficit vector to deviate from the rake specified in FaultInfo. This requires Greens functions to be computed in rake and rake-perpendicular directions. Set to ‘false’ if you want to require slip deficit in direction of rake in Faultinfo. 
(c) coseismic_slip (meters) is used to define the default recurrence time for faults. Recurrence time is coseismic_slip*1000/PrefSlipRate (PrefSlipRate in mm/yr, specified in FaultInfo.
(d) t_eq_frac is the default time since the last earthquake, normalized by recurrence time
(e) Slip patch discretization paremeters: nve_total is number of slip patches down dip on faults, pL is the nominal patch length along strike, and DLock_max is the bottom depth of the fault (also maximum locking depth).
(f) refine is the length (km) of refined fault patches used to construct the Greens Functions for the top row of slip patches. Although the inversion solves for slip deficit on larger slip patches, the refined, smaller slip patches with imposed linear slip gradients are used to construct the Greens Functions of the shallow patches to avoid strain rate singularities due to discontinuous slip gradients. (refine=1 km is recommended).
(g) Nterms is the number of terms in the Hankel transform for the viscoelastic calculations. Nterms=100 is fine for model domains of a few 100 km, but Nterms=200 shoud be used for larger model domains (about 1000 km)
(h) Viscosity model parameters: H1 is elastic thickness (km) which probably should be same as DLock_max. H2 is thickness of ‘mantle lid’ (typically higher viscosity than mantle). tR1 is relaxation time of mantle lid, and tR2 is relaxation time of mantle below depth H2.    
 
4. Various parameters/options for the inversion are specified in the INPUT SECTION of invert_strainrate_for_backslip.m:
(a) regular_grid is used only for plotting purposes. If regular_grid = true, plotting assumes a regular grid spacing of strain rate observations (in lat and long directions) with spacing of dl degrees. 
(b) If variable_rake = true, inversion will solve for slip deficit components in both the specified rake and rake-perpendicular directions (rake specified in FaultInfo) – this requires that you computed variable rake Greens functions. If variable_rake = false, inversion will fix slip deficit in direction of specified rake (this is a less expensive inversion). 
(c) rake_weight is a weighting factor placed on minimizing the amount of rake-perpendicular slip deficit rate. This allows you to control to some extent the degree to which slip deficit direction is allowed to deviate from the specified rake direction. 
(d) Set use_elastic = true to use elastic half-space Greens functions only and ignore time-dependent deformation due to viscoelastic cycle.  Set use_elastic = false to include viscoelastic earthquake cycle deformation (for faults with FaultInfo.RecurrenceTime values not specified as ‘nan’). 
(e) Set use_Gaussian_prior = true to use preferred slip rates (mean) and slip rate standard deviations in FaultInfo as gaussian priors o slip deficit rate. If bounds are included, this is a truncated Gaussian prior.
(f) If use_bounds_for_std = true, the bounds in FaultInfo will be used to construct the standard deviation for the Gaussian slip deficit rate prior. In this case, std_scale is used to compute the standard deviations as (ub-lb)*std_scale
(g) use_upper_bounds and use_lower_bounds specifiy (true/false) whether bounds are included. If use_lower_bounds = false, lower bounds on slip deficit rate are set to zero (to allow fully creeping patches). 
(h) scale_upper_bounds and scale_lower_bounds are scalar multipliers to increase or decrease bounds specified in FaultInfo (uniformly).
(i) weight_smoothing is weight placed on first-derivative smoothing of slip deficit rate across neighboring slip patches


5. Run plot_best_fitting.m to view inversion results. The INPUT SECTION you to turn on/off a number of plotting options. 

6. Run plot_principal_directions.m to view observed and model principal strain rate directions. No inputs required. 

7. Run sample_posterior.m to compute posterior distribution of slip deficit rate.  The number of samples of slip deficit rate for each slip patch can be varied (10^4 seems to work well) as well as the correlation distance. Slip deficit rate correlations for patches farther than the correlation distance will be neglected (to speed up the sampling). 

8. (optional – requires step 7) Run plot_slip_rates_along_faults.m to visualize depth-averaged slip deficit rates along specified faults. In the INPUT SECTION, provide a list of fault names (or partial fault names). Any fault name containing the name in the list will be plotted.

9. (optional – requires step 7) Run plot_moment_sources.m to visualize principal directions of moment rate sources and confidence regions. In the INPUT SECTION you can change the size of the principal direction vectors and change the legend scale bar. 

Reference:
Johnson, K. M., Wallace, L. M., Maurer, J.,Hamling, I., Williams, C., Rollins, C., et al.(2024). Inverting geodetic strain rates forslip deficit rate in complex deforming zones: An application to the New Zealandplate boundary. Journal of GeophysicalResearch: Solid Earth, 129,e2023JB027565.  https://doi.org/10.1029/2023JB027565
