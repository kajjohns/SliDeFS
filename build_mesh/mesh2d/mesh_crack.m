

%exterior nodes

node = [-10 -10;10 -10; 10 10; -10 10];
edge = [1 2; 2 3; 3 4; 4 1];


opts.kind = 'delaunay';
hfun = 4; 

% [vert,etri, tria,tnum] = refine2(node,edge,[]  ,opts,hfun) ;
% 
%  figure;
%     patch('faces',tria(:,1:3),'vertices',vert, ...
%         'facecolor','w', ...
%         'edgecolor',[.2,.2,.2]) ;
%     hold on; %axis image off;
%     patch('faces',edge(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;
%     
%     
    
    
    
%now add crack
segs = [0 -5 0 0; 0 2 2 5];
segs = [0 -5 0 0; 1 -1 2 5];
segs = [0 -5 0 0; 1 -1 2 5; 3 1 3 6];

L = 1;  %nominal patch length
length_segs = sqrt( (segs(:,1)-segs(:,3)).^2  + (segs(:,2)-segs(:,4)).^2 ); 
N = round(length_segs/L);
Px = [];
Py = [];
edge_crack = [];
for k=1:length(N)
    
    alpha = atan2(segs(k,4)-segs(k,2),segs(k,3)-segs(k,1));
    spacing=logspace(0.1,1,round(N(k)/2));
    spacing = [spacing fliplr(spacing)]
    spacing = spacing*(length_segs(k)/sum(spacing));
    spacing = [0 cumsum(spacing)];
    
    px{k} = segs(k,1) + spacing*cos(alpha);
    py{k} = segs(k,2) + spacing*sin(alpha);

    Px = [Px;px{k}'];
    Py = [Py;py{k}'];

    n = max(max([edge;edge_crack]));  %number of boundary edges
    edge_crack = [edge_crack; n+(1:length(px{k})-1)' n+(2:length(px{k}))'];

end


node_crack = [Px Py];
nodes = [node;node_crack];
edge = [edge;edge_crack]; 



part{1} = 1:size(node,1) ;
[vert,etri, tria,tnum] = refine2(nodes,[edge;edge],part  ,opts,hfun) ;

[vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;

 figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on;% axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
    
    for k=1:size(edge_crack,1)
        plot(nodes(edge_crack(k,:),1),nodes(edge_crack(k,:),2),'r')
    end
    axis equal
    
   