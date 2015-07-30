%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for lambda_ind = 1:length(lambda_range)
    for calT_ind = 1:length(calT_range)        
	for RAD_ind = 1:length(radius_range)
	    lambda = lambda_range(lambda_ind);
	    calT = calT_range(calT_ind);
	    radius = radius_range(RAD_ind);

	    savefname = sprintf('%s/solpairacc_lambda%.3f_radius%.3f_masses_calibratedT%.3f_overlap%d_o2p.mat', lossdir, lambda, radius, calT, overlap)
	    fprintf(1, 'lambda = %.2f,  calT = %.2f\n', lambda, calT);
	    if exist(savefname)
		fprintf('loading solpairacc ...');
		load(savefname);
		fprintf(' ... took %.2f s\n', toc);
		all_solpairaccs{lambda_ind}{calT_ind}{RAD_ind} = solpairacc;
	    else
		fprintf('creating solpairacc...'); tic;
		if radius ~= 1
		    all_masses = {};
		    for pim = 1:nim
			load(sprintf('%s/lambda%.3f_radius%.3f_%d_calibratedT%.3f_overlap%d_o2p.mat', massdir, lambda, radius, pim, calT, overlap));
			all_masses{pim} = masses(1:nummodes);
		    end
		end
		marginalfname = sprintf('%s/all_marginals_lambda%.3f_radius%.3f_calibratedT%.3f_overlap%d_o2p.mat', marginaldir, lambda, radius, calT, overlap);
		if exist(marginalfname)
		    load(marginalfname);
		else
		    fprintf('Creating all_node_marginals...');
		    if isempty(gcp('nocreate'))
			parpool(NUM_WORKERS);
		    end
		    node_marginals = {};
		    parfor pim = 1:nim
			if radius == 1
			    M = load(sprintf('%s/lambda0.050_radius1.000_%d_calibratedT%.3f_overlap%d_o2p.mat', marginaldir, pim, calT, overlap));
			else
			    M = load(sprintf('%s/lambda%.3f_radius%.3f_%d_calibratedT%.3f_overlap%d_o2p.mat', marginaldir, lambda, radius, pim, calT, overlap));
			end
			for solind = 1:nummodes
			    if radius ~= 1
				node_marginals{pim, solind} = M.marginals_node_beliefs{solind}(:, 2);
			    else
				node_marginals{pim, solind} = M.marginals_node_beliefs{1}(:, 2);
			    end
			end
		    end
		    save(marginalfname, 'node_marginals');
		end
		fprintf('Creating solpairacc...\n');
		%keyboard;
		if isempty(gcp('nocreate'))
		    parpool(NUM_WORKERS);
		end
		solpairacc = {};
                for pim = 1:nim
                    solpairacc{pim} = zeros(nummodes, nummodes);
                end
		parfor pim = 1:nim
		    fprintf('%d of %d\n', pim, nim);
		    [spids, N] = get_unique(double(splabels{pim}));
		    spids = spids(1:N);

		    % compute m x m loss matrix
		    acc = zeros(nummodes, nummodes);
		    Qseg = {};
		    norm_masses = zeros(nummodes, 1);
		    
		    for ps = 1:nummodes
			Qseg{ps} = {};
			norm_masses(ps) = exp(all_masses{pim}(ps) - max(all_masses{pim}));
		    end
		    norm_masses = norm_masses ./ sum(norm_masses);

		    for ps = 1:nummodes
			nnodes = length(node_marginals{pim,1})/21;
			for k = 1:21
			    if radius == 1
				Qseg{ps}{k} = label_to_seg(node_marginals{pim, 1}(21 * (0:nnodes-1) + k), double(splabels{pim}), spids);
			    else
				Qseg{ps}{k} = label_to_seg(node_marginals{pim, ps}(21 * (0:nnodes-1) + k), double(splabels{pim}), spids);
			    end
			end
		    end
		    for ps1 = 1:nummodes
			for ps2 = 1:nummodes 
			    if radius == 1
				acc(ps1, ps2) = ExpectedScore(label_to_seg(all_sols{lambda_ind}{pim}(ps1,:), double(splabels{pim}), spids), Qseg{ps2}, 21);
			    else
				acc(ps1, ps2) = ExpectedScore(label_to_seg(all_sols{lambda_ind}{pim}(ps1,:), double(splabels{pim}), spids), Qseg{ps2}, 21) * norm_masses(ps2);
			    end
			end
		    end
		    solpairacc{pim} = acc;
		end
		all_solpairaccs{lambda_ind}{calT_ind}{RAD_ind} = solpairacc;
		save(savefname,'solpairacc');
	    end
	    fprintf('...done. %.2fs\n', toc);
       end 
    end
end



