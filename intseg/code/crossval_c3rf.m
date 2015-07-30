%% cross val best settings.

ALLINDS = 1:numimages;

ntest = length(ALLINDS);
NUM_FOLDS = ntest;
SPLIT_INDS = [];
for i = 1:NUM_FOLDS
    SPLIT_INDS = [SPLIT_INDS; ((i-1)*ntest/NUM_FOLDS+1):(i*ntest/NUM_FOLDS)];
end

BEST_CURVES = zeros(NUM_FOLDS, nummodes);

for held_out_ind = 1:NUM_FOLDS
    sets_to_train_on = setdiff(1:NUM_FOLDS, held_out_ind);
    
    TRAIN_SET_INDS = SPLIT_INDS(sets_to_train_on, :);
    TRAIN_SET_INDS = TRAIN_SET_INDS(:);
    
    best_iou = -1;

    for lambda_ind = 1:length(lambda_range)
        for T_ind = 1:length(calT_range)
	    for RAD_PERCENT_ind = 1:length(radius_range)
                lambda = lambda_range(lambda_ind);  T = calT_range(T_ind);  RAD_PERCENT = radius_range(RAD_PERCENT_ind); 
  
		mbr = [];  mbrind = [];
		for pf = 1:length(TRAIN_SET_INDS)                 
		    solpairacc    = all_solpairaccs{lambda_ind, T_ind, RAD_PERCENT_ind};
		    [mbr(pf), mbrind(pf)] = max(sum(solpairacc{TRAIN_SET_INDS(pf)}, 2));
		end
		TRAIN_IOUS = all_ious{lambda_ind}(:, TRAIN_SET_INDS); 
		mbr_iou = mean(TRAIN_IOUS(sub2ind([nummodes length(TRAIN_SET_INDS)], mbrind, 1:length(TRAIN_SET_INDS))));
		if mbr_iou > best_iou
		    best_lambda_ind      = lambda_ind;
		    best_T_ind           = T_ind;
		    best_RAD_PERCENT_ind = RAD_PERCENT_ind;
		    best_iou             = mbr_iou;
		end
            end
    	end
    end

    fprintf('FOLD %d :: best params: lambda = %.2f,   T = %.2f,   R = %.2f\n', held_out_ind, lambda_range(best_lambda_ind), calT_range(best_T_ind), radius_range(best_RAD_PERCENT_ind));

    TEST_SET_INDS = SPLIT_INDS(held_out_ind, :);
    for ps = 1:nummodes
        mbr = [];  mbrind = [];
	for pf = 1:length(TEST_SET_INDS)
	    solpairacc    = all_solpairaccs{best_lambda_ind, best_T_ind, best_RAD_PERCENT_ind};
	    [mbr(pf), mbrind(pf)] = max(sum(solpairacc{TEST_SET_INDS(pf)}(1:ps, 1:ps), 2));
	end 
        TEST_IOUS = all_ious{best_lambda_ind}(:, TEST_SET_INDS);
	mbr_iou(ps) = mean(TEST_IOUS(sub2ind([nummodes length(TEST_SET_INDS)], mbrind, 1:length(TEST_SET_INDS))));
    end
    
    % add these to set of curves:
    BEST_CURVES(held_out_ind, :) = mbr_iou;
end

disp(mean(BEST_CURVES));
