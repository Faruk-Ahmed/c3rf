
fprintf('Loading image data:\n');
all_data_terms = {};   all_labels = {};   all_sparse_terms = {};     all_allL = {};
for pf = 1:numimages
    fname = flist(pf).name(1:end-4);
    fprintf('  image %d\n', pf);
    load(sprintf('%s/%s.mat', datadir, fname));
    all_data_terms{pf} = data_term;  all_labels{pf} = labels;  all_sparse_terms{pf} = sparse_term;
end

%% Compute masses:
fprintf('Computing/loading masses and marginals:\n');
MASSES = {}; MARGINALS = {};
for lambda_ind = 1:length(lambda_range)
    for T_ind = 1:length(calT_range)
	for RAD_PERCENT_ind = 1:length(radius_range) 

            lambda = lambda_range(lambda_ind);  T = calT_range(T_ind);  RAD_PERCENT = radius_range(RAD_PERCENT_ind);
            fprintf('    [lambda = %.3f,  T = %.3f,  R = %.3f]\n', lambda, T, RAD_PERCENT);

	    masses_savefname = sprintf('%s/all_masses_lambda%.2f_radius%.4f_overlap%d_T%.4f.mat', savedir_masses, lambda, RAD_PERCENT, 1-noOverlap, T);
	    marginals_savefname = sprintf('%s/marginal_data_lambda%.2f_radius%.4f_overlap%d_T%.4f.mat', savedir_marginals, lambda, RAD_PERCENT, 1-noOverlap, T);

	    if ~exist(masses_savefname, 'file') || ~exist(marginals_savefname, 'file')
		if isempty(gcp('nocreate'))
		    parpool(50);   % Set NUM_WORKERS.
		end

		for pf = 1:numimages 
		    fprintf(sprintf('        image %d\n', pf)); 
		    fname = flist(pf).name; fname = fname(1:length(fname)-4);
		    
		    data_term   = all_data_terms{pf};
		    labels      = all_labels{pf};
		    sparse_term = all_sparse_terms{pf};
		    
		    % node potentials: 
		    nnodes        = size(data_term, 2);
                    nlabels       = size(data_term, 1);
		    ne            = data_term/T;
		    ne([1, 2], :) = ne([2, 1], :);

		    %% edge potentials:
		    [node1, node2, wt] = find(triu(sparse_term));
		    wt                 = wt/T;
		    nedges             = length(wt);
		    el                 = [node1 node2]';
		    ee                 = zeros(4,nedges);
		    ee(2,:)            = wt;
		    ee(3,:)            = wt;
		    
		    edgePot            = zeros(2, 2, nedges);
		    for e = 1:nedges
			edgePot(:, :, e) = [0 -wt(e); -wt(e) 0];
		    end

		    %% other stuff:
		    nStates          = 2*ones(1, nnodes);
		    adjacency_matrix = (sparse_term > 0);
		    nodePot          = -ne';
		    
		    %% settings for UGM inference:
		    max_iter   = 1000;  
		    useMex     = 0;
		    edgeStruct = UGM_makeEdgeStruct(el', nStates, nnodes, useMex, max_iter);

                    %% Load divsols:
		    divsol_savefname = sprintf('%s/%s_%.2f.mat', divsoldir, fname, lambda);
		    divsolseg_savefname = sprintf('%s/%s_segs_mmodes_bound_lambda%.2f.mat', divsoldir, fname, lambda);
                    if exist(divsol_savefname, 'file') && exist(divsolseg_savefname, 'file')
                        fprintf('          --  loading divsols\n');
                        DIVSOLS = load(divsol_savefname);
                        all_allL{pf} = DIVSOLS.allL;
                    else
                        gt = imread(sprintf('%s/%s.png', fullfile(datadir, 'gtdir'), fname));
                        allL = zeros(nummodes, nnodes);
                        sol_iou = zeros(nummodes, 1);
                        seg  = {};
                        fprintf('          --  divsols not present. Computing them.\n');
                        divne = ne * T;
                        
                        for ps = 1:nummodes
                            [L en lb]   = perform_inference(divne, el, ee * T, 'trw');

                            allL(ps,:)  = L';
                            seg{ps}     = label2seg(L, labels);
                            [~, ~, ~, sol_iou(ps), ~] = computeStats(seg{ps}, gt);
                            
                            bedges = find(L(el(1,:)) ~= L(el(2,:)));
	                    bnodes = unique(el(:,bedges));
	                    inds = sub2ind([nlabels nnodes],L(bnodes)+1,bnodes);

                            divne(inds) = divne(inds) + lambda;
                        end
                        all_allL{pf} = allL;
                        save(divsol_savefname, 'allL', 'sol_iou');
                        save(divsolseg_savefname, 'seg');
                    end

		    %%
		    sols        = 1:nummodes;
		    masses      = zeros(1, length(sols));
		    n_marginals = {};
		    e_marginals = {};

		    parfor s_ind = 1:length(sols)
			solind = sols(s_ind);
			SOLUTION       = all_allL{pf}(solind, :);
			HAMMING_RADIUS = floor(RAD_PERCENT*nnodes);
			if noOverlap
			    intersol_hamming_dist = [];
			    for othersolind = setdiff(1:length(sols), s_ind)
				intersol_hamming_dist(end+1) = sum(all_allL{pf}(sols(s_ind),:) ~= all_allL{pf}(sols(othersolind),:));
			    end
			    HAMMING_RADIUS = min(HAMMING_RADIUS, floor(0.5*min(intersol_hamming_dist)));
			end
                        fprintf('        * Computing masses and marginals for image %d, divsol %d\n', pf, solind);			
			[MASS, node_marginals, edge_marginals]  = get_mass(nodePot, edgePot, edgeStruct, SOLUTION, HAMMING_RADIUS, 0, pf, lambda, solind);
				
			masses(s_ind)      = MASS;
			n_marginals{s_ind} = node_marginals;
			e_marginals{s_ind} = edge_marginals;
		    end
		    %%
		    all_masses{pf}         = masses;
		    all_node_marginals{pf} = n_marginals;
		end
		fprintf('    Saving masses and marginals...\n');
		save(masses_savefname, 'all_masses');
		save(marginals_savefname, 'all_node_marginals');
		%MASSES{lambda_ind, T_ind, RAD_PERCENT_ind}.all_masses = all_masses;
		%MARGINALS{lambda_ind, T_ind, RAD_PERCENT_ind}.all_node_marginals = all_node_marginals;
            end
	end
    end
end
all_allL = {};
all_ious = {};
for lambda_ind = 1:numel(lambda_range)
    all_allL{lambda_ind} = {};
    all_ious{lambda_ind} = zeros(nummodes, numimages);
    for pf = 1:numimages
	fname               = flist(pf).name(1:end-4);
	divsol_savefname    = sprintf('%s/%s_%.2f.mat', divsoldir, fname, lambda);
	divsolseg_savefname = sprintf('%s/%s_segs_mmodes_bound_lambda%.2f.mat', divsoldir, fname, lambda);

	DIVSOLS                     = load(divsol_savefname);
	all_allL{lambda_ind}{pf}    = DIVSOLS.allL;
	all_ious{lambda_ind}(:, pf) = DIVSOLS.sol_iou;
    end
end
