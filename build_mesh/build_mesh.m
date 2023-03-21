
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT SECTION

%mat file for saving triangular mesh
savename = 'CSAF_mesh';


%nominal slip patch length for fault discretization
%(note) this is required to make the mesh so that triangular nodes
%correspond with endpoints of fault patches
%This will need to be the same as the nominal patch length specified in
%build_backslip_GreensFunctions.m
pL = 15; 

%largest triangle side length
hmax = 100;

% load fault info (only need SegEnds for building mesh
% this needs to be a mat file that contains a cell array of FaultInfo
% structures for each fault in the model. For the k-th fault, FaultInfo{k}
% contains
%     FaultInfo{k}.SegEnds,  (Nx4 matrix of N fault segment coordinates (Cartesian) of
%       format [x1 x2 y1 y2] where (x1,y1) is first endpoint and (x2,y2)
%       is the other endpoint

load fault_info_CSAF.mat


%edges of model domain (specify x and y limits of model domain)
model_xmin = -500; model_xmax = -200;
model_ymin = -100; model_ymax = 150;

%% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath mesh2d
addpath tools


%put all segments together into a single matrix
SegEnds = [];
for k=1:length(FaultInfo)
    SegEnds = [SegEnds; FaultInfo{k}.SegEnds];
end

%build "edges" of fault patch model along fault segments for use in constructin the triangular mesh
make_patches_mesh

%need to remove repeated nodes because
%meshing algorithm doesn't like repeats
node = round(node,4);
[C,ia,ic] = unique(node,'rows','stable');
node = C;


%need to renumber nodes values in edges to be consistent 
%with unique node values
new_edge = edge;
for k=1:length(ic)
    new_edge(edge==k)=ic(k);
end

edge = new_edge;

%edge = edge(1:350,:);

node = node(1:max(max(edge)),:);

%boundary nodes
x1 = model_xmin; x2 =model_xmax;
y1 = model_ymin; y2 = model_ymax;

bn = [x1 y1;x1 y2;x2 y2;x2 y1];
node = [bn; node];

edge = [1 2;2 3;3 4;4 1;edge+4];


 part{1} = [1,2,3,4] ;
    
%---------------------------------------------- do size-fun.
   
 
%---------------------------------------------- do mesh-gen.
    opts.kind = 'delaunay' ;
   
   [vert,etri, tria,tnum] = refine2(node,edge,part,opts,hmax) ;

    [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
  
   
   
   figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; 
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',2.0) ;

    axis equal


%     
%     %plot gps sations
%     
%     A = load('gnss_data.txt');
%     origin = [mean(A(:,2)) mean(A(:,1))];
%     xy_gps = llh2local(A(:,1:2)',fliplr(origin))';
%     plot(xy_gps(:,1),xy_gps(:,2),'r.')
%     
%     
    tri=tria;
    nodes = vert;
    
    
    plot(node(:,1),node(:,2),'ro')
    
    title('triangular mesh (red circles are slip patch ends)')
    
    
    patch_stuff=make_triangular_patch_stuff(tri,[nodes 0*nodes(:,1)]);
    
    tri_areas = patch_stuff.area_faces;
    tri_centroids = patch_stuff.centroids_faces;
    
 save(savename, 'tri', 'nodes', 'tri_areas', 'tri_centroids')    
    