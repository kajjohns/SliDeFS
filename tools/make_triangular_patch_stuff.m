function patch_stuff=make_triangular_patch_stuff(tri,p)

strikevec_faces=[];
strike_faces=[];
dipvec_faces=[];
dip_faces=[];
strike_faces=[];
dip_faces=[];
centroids_faces=[];
normal_faces=[];
area_faces=[];


for j=1:size(tri,1)

    temp1{1}=p;
    temp2{1}=tri(j,:);

    %get strike and dip and normal vectors and area for each triangular patch
    vec1=p(temp2{1}(1),:)-p(temp2{1}(2),:);
    vec2=p(temp2{1}(3),:)-p(temp2{1}(2),:);
    cross_face=cross(vec1,vec2);
    veclength=sqrt(cross_face(1).^2+cross_face(2).^2+cross_face(3).^2);
    normal=cross_face/veclength;
    strikevec=[1 -normal(1)/normal(2) 0]/norm([1 -normal(1)/normal(2) 0]);
    dipvec=cross(normal,strikevec);

    %change signs
    if dipvec(3)>0
        dipvec=-dipvec;
    end
    if normal(3)<0
        normal=-normal;
    end
    strikevec=cross(normal,dipvec);



    normal_faces=[normal_faces;normal];
    strikevec_faces=[strikevec_faces;strikevec];
    dipvec_faces=[dipvec_faces;dipvec];
    strike_faces=[strike_faces;90-atan2(strikevec(2),strikevec(1))*180/pi];
    dip_faces=[dip_faces;abs(atan(dipvec(3)/sqrt(dipvec(1)^2+dipvec(2)^2))*180/pi)];
    centroids_faces=[centroids_faces;[mean(temp1{1}(temp2{1},1)) mean(temp1{1}(temp2{1},2)) mean(temp1{1}(temp2{1},3))]];
    area_faces=[area_faces;.5*abs(norm(cross(vec1,vec2)))];


end %j


patch_stuff.strikevec_faces=strikevec_faces;
patch_stuff.strike_faces=strike_faces;
patch_stuff.dipvec_faces=dipvec_faces;
patch_stuff.dip_faces=dip_faces;
patch_stuff.strike_faces=strike_faces;
patch_stuff.dip_faces=dip_faces;
patch_stuff.centroids_faces=centroids_faces;
patch_stuff.normal_faces=normal_faces;
patch_stuff.area_faces=area_faces;
