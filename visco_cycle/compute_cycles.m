
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

        
        Ld = lower_locking_depths(k)*ones(numsegs,1);
        ULd = upper_locking_depths(k)*ones(numsegs,1);
        Ts = T(k)*ones(numsegs,1);
        teqs = teq(k)*ones(numsegs,1);
        dips = Dips(k)*ones(numsegs,1);
        
        
        [Gvel,Gstrain] = Get_interseismic_strainrates_cycle(Segs,dips,Ld,ULd,...
            xystats,xs,ys,xy_visco,H1,H2,tR1,tR2,Ts,teqs,get_SS,get_DS,Nterms);

        Ve_ss(:,k) = sum(Gvel.ss_ve,2);  %sum all segments
        Vn_ss(:,k) = sum(Gvel.ss_vn,2);
        Vu_ss(:,k) = sum(Gvel.ss_vu,2);

        Ve_ds(:,k) = sum(Gvel.ds_ve,2);
        Vn_ds(:,k) = sum(Gvel.ds_vn,2);
        Vu_ds(:,k) = sum(Gvel.ds_vu,2);

        Ve_ss_elastic(:,k) = sum(Gvel.ss_ve_e,2);
        Vn_ss_elastic(:,k) = sum(Gvel.ss_vn_e,2);
        Vu_ss_elastic(:,k) = sum(Gvel.ss_vu_e,2);

        Ve_ds_elastic(:,k) = sum(Gvel.ds_ve_e,2);
        Vn_ds_elastic(:,k) = sum(Gvel.ds_vn_e,2);
        Vu_ds_elastic(:,k) = sum(Gvel.ds_vu_e,2);

        Exx_ss(:,k) = sum(Gstrain.ss_Exx,2);
        Exy_ss(:,k) = sum(Gstrain.ss_Exy,2);
        Eyy_ss(:,k) = sum(Gstrain.ss_Eyy,2);

        Exx_ds(:,k) = sum(Gstrain.ds_Exx,2);
        Exy_ds(:,k) = sum(Gstrain.ds_Exy,2);
        Eyy_ds(:,k) = sum(Gstrain.ds_Eyy,2);

        Exx_ss_elastic(:,k) = sum(Gstrain.ss_e_Exx,2);
        Exy_ss_elastic(:,k) = sum(Gstrain.ss_e_Exy,2);
        Eyy_ss_elastic(:,k) = sum(Gstrain.ss_e_Eyy,2);

        Exx_ds_elastic(:,k) = sum(Gstrain.ds_e_Exx,2);
        Exy_ds_elastic(:,k) = sum(Gstrain.ds_e_Exy,2);
        Eyy_ds_elastic(:,k) = sum(Gstrain.ds_e_Eyy,2);
        