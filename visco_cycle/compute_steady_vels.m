
        Segs = SegEnds_xy{k};
        numsegs = size(Segs,1);

        if Rakes(k)==0 | Rakes(k)==180
            get_DS = false(numsegs,1);
        else
            get_DS = true(numsegs,1);
        end

        if Rakes(k)==90 | Rakes(k)==-90
            get_SS = false(numsegs,1);
        else
            get_SS = true(numsegs,1);
        end

        dips = Dips(k)*ones(numsegs,1);
        
       
        %long-term (steady-state) velocities
        [Gss, Gds] = Get_steady_vels(Segs,dips,xystats,xy_visco,H1,get_SS,get_DS);

        Ve_steady_ss(:,k) = sum(Gss.east,2);  %sum all segments
        Vn_steady_ss(:,k) = sum(Gss.north,2);
        Vu_steady_ss(:,k) = sum(Gss.up,2);

        Ve_steady_ds(:,k) = sum(Gds.east,2);
        Vn_steady_ds(:,k) = sum(Gds.north,2);
        Vu_steady_ds(:,k) = sum(Gds.up,2);

        