VOCinit;
nsol = nummodes;

if ~exist(massdir, 'dir')
    unix(sprintf('mkdir %s', save_folder_name));
end
if ~exist(marginaldir, 'dir')
    unix(sprintf('mkdir %s', marginal_save_folder_name));
end

flistfile = sprintf(VOCopts.seg.imgsetpath,VOCopts.testset);

if exist(sprintf('%s/%s_stats.mat',savedir,VOCopts.dataset),'file')
    load(sprintf('%s/%s_stats.mat',savedir,VOCopts.dataset));
else
    for pts = 1:length(testsets)
        testset = testsets{pts}

        imnames = textread(sprintf(VOCopts.seg.imgsetpath,testset),'%s\n');
        imnums = char(imnames); % convert to char 2D array first
        imnums = str2num([imnums(:,1:4) imnums(:,6:end)]); % drop the _ and convert to number

        nim = length(imnames);
        sizes = zeros(nim,2);

        for pim=1:length(imnames)
            if mod(pim,50) == 0
                fprintf('%d...\t',pim);
            end

            im = imread(sprintf(VOCopts.imgpath,imnames{pim}));
            sizes(pim,:) = [size(im,1) size(im,2)];
        end
        imsizes{pts} = sizes;
    end
    save(sprintf('%s/%s_stats.mat',savedir,VOCopts.dataset),'imsizes');
end

imnames = textread(flistfile,'%s\n');
imnums = char(imnames);
imnums = str2num([imnums(:,1:4) imnums(:,6:end)]);

if exist(sprintf('%s/%s%s_gtsegsen.mat',savedir,VOCopts.dataset,VOCopts.testset),'file')
    load(sprintf('%s/%s%s_gtsegsen.mat',savedir,VOCopts.dataset,VOCopts.testset));
    for pim = 1:nim
        gtsegs{pim} = double(gtsegs{pim});
    end
else
    for pim=1:length(imnames)
        gtsegs{pim} = double(imread(fullfile(clsimgpath, sprintf('%s.png',imnames{pim}))));
    end
    save(sprintf('%s/%s%s_gtsegsen.mat',savedir,VOCopts.dataset,VOCopts.testset),'gtsegs');
end
gtsegs = gtsegs(1:nim);
sol_ens    = {};
all_masses = {};

%% load data:
fprintf('Loading image data:\n');
for pim = 1:nim
    if mod(pim, 50) == 0
        fprintf('  image %d\n', pim);
    end
    D = load(sprintf('%s/data_%d.mat', divsoldir, pim));
    splabels{pim} = double(D.sp);
    crfnes{pim} = D.ne;
    crfels{pim} = D.el;
    crfnedges{pim} = D.nedges;
    crfpottsterms{pim} = D.pottsterm;
end

fprintf('Loading divsols:\n');
for lambda_ind = 1:numel(lambda_range)
    lambda = lambda_range(lambda_ind);
    fprintf('  lambda = %.4f:\n', lambda);

    load(sprintf('%s/all_o2p_divsols_lambda%.3f.mat', divsoldir, lambda));
    remove_dupes = 1;
    if remove_dupes
        fprintf('    Creating sol set without duplicates\n');
        dupes_removed_sols = {};
        for pim = 1:nim
            dupes_removed_sols{pim} = unique(all_sols{pim}(1:50,:), 'rows', 'stable');
            if size(dupes_removed_sols{pim}, 1) < nummodes
                dupes_removed_sols{pim} = repmat(dupes_removed_sols{pim}, nummodes, 1);
            end
            dupes_removed_sols{pim} = dupes_removed_sols{pim}(1:nummodes, :);
        end
        for pim = 1:nim
            all_sols_new{lambda_ind}{pim} =  double(dupes_removed_sols{pim});
        end
    end
    all_sols = all_sols_new;
    
    savefname_solsegs = sprintf('%s/solsegs_nodupes_%.4f_M%d.mat', divsoldir, lambda, nummodes);
    solsegs = {};
    if ~exist(savefname_solsegs, 'file')
        fprintf('    Creating solsegs:\n');
	for pim = 1:nim
            fprintf('      image %d \n', pim);
            [spids, NN] = get_unique(double(splabels{pim}));
            spids = spids(1:NN);
            for ps = 1:nummodes
                  solsegs{pim, ps} = double(label_to_seg(all_sols_new{lambda_ind}{pim}(ps, :), double(splabels{pim}), spids));
           end
       end
       fprintf('Saving solsegs...\n');
       save(savefname_solsegs, 'solsegs', '-v7.3');
       fprintf('\n');
   else
      fprintf('  Loading solsegs ...\n');
      load(savefname_solsegs);
   end

    all_solsegs{lambda_ind} = solsegs;

    for RAD_PERCENT = radius_range
        parfor pim = 1:nim
	    compute_masses(lambda, RAD_PERCENT, noOverlap, calT_range, nummodes, pim, all_sols{lambda_ind}, massdir, marginaldir, divsoldir);
        end
    end
end
