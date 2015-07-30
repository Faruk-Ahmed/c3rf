
SET_INDS = 1:numimages;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntest = length(SET_INDS);

fprintf('Computing/loading pairwise losses:\n');
for lambda = lambda_range
    for T = calT_range
        for RAD_PERCENT = radius_range

	    fprintf('    [lambda = %.3f,  T = %.3f,  R = %.3f]\n', lambda, T, RAD_PERCENT);
	    %% compute pairwise loss
	    savefname = sprintf('%s/solpairacc_mmodes_bound_lambda%.2f_T%.2f_radius%.4f_overlap%d_c3rf.mat', savedir_pairwise_losses, lambda, T, RAD_PERCENT, 1-noOverlap);

	    if exist(savefname, 'file')
		load(savefname);
	    else
                solpairacc = {};
                for pf = 1:numimages
                    solpairacc{pf} = zeros(nummodes, nummodes);
                end
		%% load masses and marginals:
		MARGINALS = load(sprintf('%s/marginal_data_lambda%.2f_radius%.4f_overlap%d_T%.4f.mat', savedir_marginals, lambda, RAD_PERCENT, 1-noOverlap, T));
		MASSES = load(sprintf('%s/all_masses_lambda%.2f_radius%.4f_overlap%d_T%.4f.mat', savedir_masses, lambda, RAD_PERCENT, 1-noOverlap, T));

		if isempty(gcp('nocreate'))
		    parpool(50);
		end
		parfor pf = 1:numimages
		    fprintf(1, '      * image %d \n', pf);
		    fname = flist(pf).name; fname = fname(1:length(fname)-4);
		    DIVSOLS = load(sprintf('%s/%s_%.2f.mat', divsoldir, fname, lambda));
		    IMGDATA = load(sprintf('%s/%s.mat', datadir, fname));
		    SEG     = load(sprintf('%s/%s_segs_mmodes_bound_lambda%.2f.mat',divsoldir, fname, lambda));

		    % compute m x m loss matrix
		    iou = zeros(nummodes, nummodes);		    
		    margSegs = {};
		    for ps = 1:nummodes
			if RAD_PERCENT == 1
			    Q = MARGINALS.all_node_marginals{pf}{1}(:, 2);
			else
			    Q = MARGINALS.all_node_marginals{pf}{ps}(:, 2);
			end
			margSegs{ps} = label2seg(Q, IMGDATA.labels);
		    end
		    
		    for ps1 = 1:nummodes
			for ps2 = 1:nummodes
			    
			    Q1 = margSegs{ps2};
			    Q0 = 1 - Q1;
			    
			    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			    prediction0 = SEG.seg{ps1} == 0;
			    expectedIntersection0 = min(Q0, prediction0);
			    expectedUnion0        = max(Q0, prediction0);

			    prediction1 = SEG.seg{ps1} == 1;
			    expectedIntersection1 = min(Q1, prediction1) ;
			    expectedUnion1        = max(Q1, prediction1);

			    WeightedIOU = 0.5*(sum(expectedIntersection1(:))/sum(expectedUnion1(:)) + sum(expectedIntersection0(:))/sum(expectedUnion0(:)));
			    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			    iou(ps1, ps2) = WeightedIOU * exp(MASSES.all_masses{pf}(ps2) - max(MASSES.all_masses{pf}));
			end
		    end
		    solpairacc{pf} = iou;
		end
		fprintf('\n');
		save(savefname, 'solpairacc');
	    end
	    all_solpairaccs{find(lambda_range == lambda), find(calT_range == T), find(radius_range == RAD_PERCENT)} = solpairacc;
	end
    end
end

