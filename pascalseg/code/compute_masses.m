function compute_masses(lambda, RAD_PERCENT, noOverlap, Trange, nummodes, ind, all_sols, massdir, marginaldir, divsoldir)

for calT = Trange
    fprintf('  * Computing masses/marginals for image %d\n', ind);
    ALPHA = 0.3;

    allL = all_sols{ind} + 1;

    flag = 0;
    masses = [];
    marginals_node_beliefs = {};
    marginals_edge_beliefs = {};

    save_fname = sprintf('%s/lambda%.3f_radius%.3f_%d_calibratedT%.3f_overlap%d_o2p.mat', massdir, lambda, RAD_PERCENT, ind, calT, 1-noOverlap);
    save_marginals_fname = sprintf('%s/lambda%.3f_radius%.3f_%d_calibratedT%.3f_overlap%d_o2p.mat', marginaldir, lambda, RAD_PERCENT, ind, calT, 1-noOverlap);

    if ~exist(save_fname, 'file') || ~exist(save_marginals_fname, 'file')
	ne = [];    
	DATA = load(sprintf('%s/data_%d.mat', divsoldir, ind));

	nlabels = size(DATA.ne, 1);
	nnodes = size(DATA.ne, 2);
	nodePot = -DATA.ne'/calT;
	adj = zeros(nnodes, nnodes);
	for s_i = 1:DATA.nedges
	    adj(DATA.el(1,s_i), DATA.el(2,s_i)) = 1;
	end
	adj = adj + adj';

	useMex = 0;    maxIter = 1;
	nStates = nlabels;

	nodePotBin = zeros(nlabels*nnodes,2);

	for s_i = 1:nlabels
	    nodePotBin( nlabels*(0:nnodes-1)+s_i, 2) = nodePot(:,s_i);
	    nodePotBin( nlabels*(0:nnodes-1)+s_i, 1) = 0;
	end

	elBin = [];
	[p, q] = meshgrid(1:nlabels, 1:nlabels);
	pairs = [p(:) q(:)];
	multiEdgeInd = [];
	for s_i = 1:DATA.nedges
	    elBin = [elBin; ((DATA.el(1,s_i)-1)*nlabels + pairs(:,1)) ((DATA.el(2,s_i)-1)*nlabels + pairs(:,2))];
	    multiEdgeInd = [multiEdgeInd; repmat(s_i, size(pairs,1), 1)];
	end

	nedgesBin = size(elBin, 1);

	adjBin = zeros(nnodes*nlabels, nnodes*nlabels);
	for s_i = 1:nedgesBin
	    adjBin(elBin(s_i, 1), elBin(s_i, 2)) = 1;
	end

	adjBin = adjBin + adjBin';

	edgeStructBin = UGM_makeEdgeStruct(elBin, 2, 0, 1, nnodes*nlabels);

	edgePotBin = ones(2, 2, nedgesBin);
	for e = 1:nedgesBin
	    edgePotBin(:, :, e) = [0 0; 0 -1e-4*DATA.pottsterm(multiEdgeInd(e))/calT];
	end

	intersolution_hamming_distance = zeros(nummodes, nummodes);
	for sol1 = 1:nummodes
	    for sol2 = 1:nummodes
		intersolution_hamming_distance(sol1, sol2) = sum(allL(sol1, :) ~= allL(sol2, :));
	    end
	end
	intersolution_hamming_distance(find(eye(nummodes))) = Inf;

	for solind = 1:nummodes
	    
	    %% Initialize things:
	    useDamping = 0;
	    messages_old = zeros(2, 2*nedgesBin);    maximize = 0;
	    old_node_beliefs = zeros(size(nodePotBin));
	    nodeBelOld = zeros(size(nodePotBin));
	    nodeBelOld_inner = zeros(size(nodePotBin));
	    MSGS_FOR_UP_CARD_TREE_OLD = zeros(size(nodePotBin));

	    NODE_POT = nodePotBin;
	    MSG_FROM_NODE_POT = exp(my_norm_log(nodePotBin));
	    MSGS_FROM_NODE = MSG_FROM_NODE_POT;
	    MSGS_FROM_XOR = ones(size(MSGS_FROM_NODE));
	    MSGS_FROM_CARD = ones(size(MSGS_FROM_NODE));
	    iter = 0;

	    HAMMING_RADIUS = floor(RAD_PERCENT * nnodes);
	    if noOverlap
		HAMMING_RADIUS = min(HAMMING_RADIUS, floor(0.5*min(intersolution_hamming_distance(solind, :))));
	    end
	    
	    L = allL(solind,:)';
	    seg_inds = unique(L);

	    nNeighbors = zeros(size(nodePotBin,1),1);
	    nNeighbors(21*(0:nnodes-1)+L') = nNeighbors(21*(0:nnodes-1)+L') + 1;

	    while 1
		iter = iter + 1;
		calc_margs = 1;

		if iter > 20
		    useDamping = 1;
		    ALPHA = min(ALPHA + iter*0.001, 0.5);
		end
		
		
		MSGS_DOWN_TREES = (MSGS_FROM_CARD ./ repmat(sum(MSGS_FROM_CARD,2),1,2)) .* (MSGS_FROM_XOR ./ repmat(sum(MSGS_FROM_XOR,2),1,2));
		MSGS_DOWN_TREES = MSGS_DOWN_TREES ./ repmat(sum(MSGS_DOWN_TREES,2), 1, 2);
		    
		[nodeBel, edgeBel, logZ_SP, new_msg, MSGS_FOR_UP_TREE] = Infer_LBP_modified(MSGS_FROM_NODE, ...
											    edgePotBin, ...
											    edgeStructBin, ...
											    NODE_POT, ...
											    MSGS_DOWN_TREES, ...
											    nNeighbors, ...
											    calc_margs, ...
											    useDamping, ...
											    messages_old, ...
											    ALPHA, ...
											    MSG_FROM_NODE_POT);
		
		%% Impose XOR:
		MSGS_FOR_UP_XOR_TREE = MSGS_FOR_UP_TREE .* MSGS_FROM_CARD;
		MSGS_FOR_UP_XOR_TREE = MSGS_FOR_UP_XOR_TREE ./ repmat(sum(MSGS_FOR_UP_XOR_TREE, 2), 1, 2);
		[MSGS_FROM_XOR, logZ_XOR] = get_msg_from_XOR_factor(MSGS_FOR_UP_XOR_TREE, DATA.nnodes, nlabels, 1);

		if sum(abs(nodeBel(:) - nodeBelOld_inner(:))) < 1e-3
		   break;
		end
		nodeBelOld_inner = nodeBel;
		MSGS_FROM_NODE = MSG_FROM_NODE_POT .* MSGS_FROM_XOR .* MSGS_FROM_CARD;
		MSGS_FROM_NODE = exp(my_norm_log(my_log(MSGS_FROM_NODE))) ;

		MSGS_FROM_CARD = ones(size(NODE_POT));

		MSGS_FOR_UP_CARD_TREE = MSGS_FOR_UP_TREE .* MSGS_FROM_XOR;
		MSGS_FOR_UP_CARD_TREE = exp(my_norm_log(my_log(MSGS_FOR_UP_CARD_TREE)));

		LABEL_INDS = nlabels*(0:nnodes-1)+L';
		[~, messages_label, logZ_card] = sum_cardinality(int32(length(LABEL_INDS)), ...
								  my_log(MSGS_FOR_UP_CARD_TREE(LABEL_INDS, 2)), ...
								  my_log(MSGS_FOR_UP_CARD_TREE(LABEL_INDS, 1)), ...
								  int32(calc_margs), ...
								  int32(HAMMING_RADIUS));

		MSGS_FROM_CARD(LABEL_INDS,:) = exp(messages_label(:, [2,1]));
		MSGS_FROM_NODE = MSG_FROM_NODE_POT .* MSGS_FROM_XOR .* MSGS_FROM_CARD;
		MSGS_FROM_NODE = exp( my_norm_log(my_log(MSGS_FROM_NODE)) ) ;

		messages_old = new_msg;
	    end
	    if (flag == 1) || (iter > 5000)
		fprintf('Early stopping with flag = %d', flag);
		break
	    end

	    masses(solind) = logZ_card + logZ_SP + logZ_XOR;
	    marginals_node_beliefs{solind} = nodeBel;
	end
	save(save_fname, 'masses');
	save(save_marginals_fname, 'marginals_node_beliefs');
    end
end

