
%Dip_Dir (1=NE,2=SE,3=SW,4=NW,5=N,6=E,7=S,8=W)

for j=1:size(SegEnds_xy{k},1)

    if strmatch('NE',DipDir{k})
        testdir = [1 1];
    elseif strmatch('SE',DipDir{k})
        testdir = [1 -1];
    elseif strmatch('SW',DipDir{k})
        testdir = [-1 -1];
    elseif strmatch('NW',DipDir{k})
        testdir = [-1 1];
    elseif strmatch('N',DipDir{k})
        testdir = [0 1];
    elseif strmatch('E',DipDir{k})
        testdir = [1 0];
    elseif strmatch('S',DipDir{k})
        testdir = [0 -1];
    elseif strmatch('W',DipDir{k})
        testdir = [-1 0];       
    else %vertical
        testdir = [0 0];
    end
    
    %dot test direction into normal
    dt = dot(SegNorm{k}(j,:),testdir);
    
    %flip direction if dt<0 (testdir is in wrong quadrant)
    if dt<0
       SegEnds_xy{k}(j,:) = [SegEnds_xy{k}(j,3:4) SegEnds_xy{k}(j,1:2)]; %flip ends
       SegEnds_llh{k}(j,:) = [SegEnds_llh{k}(j,3:4) SegEnds_llh{k}(j,1:2)]; %flip ends
    end

end

