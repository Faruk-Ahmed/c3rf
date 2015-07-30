
NCLASSES = VOCopts.nclasses;

if 1
    NUM_PERMS = 5;
    NUM_FOLDS = 10;
    ALLINDS = 1:nim;
    SPLIT_INDS = {};
    for perm = 1:NUM_PERMS
        CROSS_INDS = crossvalind('Kfold', nim, NUM_FOLDS);
        SPLIT_INDS{perm} = {};
        for i = 1:NUM_FOLDS
            SPLIT_INDS{perm}{i} = ALLINDS(CROSS_INDS == i);
        end
    end
else
    NUM_FOLDS = 10;    NUM_PERMS = 5;
    load(sprintf('SPLIT_INDS_NUMPERMS_%d_NUMFOLDS_%d.mat', NUM_PERMS, NUM_FOLDS));
    disp('loading SPLIT_INDS');
end

SOLSEGS_TEST_INDS = {};
for ps = 1:nummodes
    SOLSEGS_TEST_INDS{ps} = {};
end

classwise_accuracies = {};
mass_accs = {};
for perm = 1:NUM_PERMS
    mass_accs{perm}  = zeros(nummodes, 1);
    classwise_accuracies{perm} = [];
end

for perm = 1:NUM_PERMS
    for held_out_ind = 1:NUM_FOLDS
        held_out_ind
	TRAIN_SET_INDS = [];
        for i = setdiff(1:NUM_FOLDS, held_out_ind);
            TRAIN_SET_INDS = [TRAIN_SET_INDS  SPLIT_INDS{perm}{i}];
        end
    
        best_mass_mc_over_MAP = 0;
        
        best_M_delta = 0; best_M_mass_mc = 0;  best_M_mass = 0;  best_M_full_marginals = 0;
        for lambda_ind = 1:numel(lambda_range)
            for calT_ind = 1:numel(calT_range)
		for RAD_ind = 1:numel(radius_range)
		    lambda = lambda_range(lambda_ind);
		    calT   = calT_range(calT_ind);
		    RAD_PERCENT = radius_range(RAD_ind);
		    [inds_masses_train] = mbr_expected_loss(TRAIN_SET_INDS, [], lambda, RAD_PERCENT, nummodes, overlap, calT, all_solsegs{lambda_ind}, all_solpairaccs{lambda_ind}{calT_ind}{RAD_ind});
		    %%%%%%%%%%%%%%%%%%%%%% evaluate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    predsegs_masses = {};   gtsegs_CV = {};
		    for pi = 1:length(TRAIN_SET_INDS)
			predsegs_masses{pi} = all_solsegs{lambda_ind}{TRAIN_SET_INDS(pi), inds_masses_train(pi)};
			gtsegs_CV{pi} = gtsegs{TRAIN_SET_INDS(pi)};
		    end

		    accuracies_train_masses = VOCevalseg_c(NCLASSES, predsegs_masses, gtsegs_CV);
		    
		    max_over_MAP_mass_mc = mean(accuracies_train_masses);
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		    if max_over_MAP_mass_mc > best_mass_mc_over_MAP
			 best_lambda_mass_ind = lambda_ind;  best_calT_mass_ind = calT_ind;
			 best_mass_mc_over_MAP = max_over_MAP_mass_mc;
			 best_RAD_ind = RAD_ind;
		    end
		end
	    end
        end

        fprintf('best_calT = %f,   best_lambda = %f,   best_RAD = %f\n', calT_range(best_calT_mass_ind), lambda_range(best_lambda_mass_ind), radius_range(best_RAD_ind));

        TEST_SET_INDS = SPLIT_INDS{perm}{held_out_ind};
    
        for ps = 1:nummodes 
            [inds_masses_test] = mbr_expected_loss(TEST_SET_INDS, ...
                            [],lambda_range(best_lambda_mass_ind),...						
                 	    radius_range(best_RAD_ind), ps, overlap,...
                 	    calT_range(best_calT_mass_ind),...
                 	    all_solsegs{best_lambda_mass_ind},...
                 	    all_solpairaccs{best_lambda_mass_ind}{best_calT_mass_ind}{best_RAD_ind});
    
            for pind = 1:length(TEST_SET_INDS)
                pim = TEST_SET_INDS(pind);
                SOLSEGS_TEST_INDS{ps}{pim} = all_solsegs{best_lambda_mass_ind}{pim, inds_masses_test(pind)};
            end
        end
    end

    for ps = 1:nummodes
        accuracies = VOCevalseg_c(NCLASSES,SOLSEGS_TEST_INDS{ps}, gtsegs);
	mass_accs{perm}(ps) = mean(accuracies);
        fprintf('%.2f\n', mean(accuracies));
	if ps == nummodes
            classwise_accuracies{perm} = accuracies;
	end
    end
end

D = cell2mat(mass_accs)
mean(D, 2)


