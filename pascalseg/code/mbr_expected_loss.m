function [inds_masses] = mbr_expected_loss(SET_INDS, Trange, lambda, radius, nummodes, overlap, calT, solsegs, solpairacc)
    ps = nummodes
    mbr_masses = [];  mbrind_masses = [];
    for pind=1:length(SET_INDS)
	pim = SET_INDS(pind);
	% because we have accuracy not loss
	[mbr_masses(pind) mbrind_masses(pind)] = max(sum(solpairacc{pim}(1:ps,1:ps)'));
    end
    inds_masses = mbrind_masses;
   
