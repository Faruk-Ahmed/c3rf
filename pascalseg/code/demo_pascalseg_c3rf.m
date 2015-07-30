clear; close; clc;
dbstop if error;

%% set paths:
datadir = '../data';
divsoldir = fullfile(datadir, 'divsols');
savedir = fullfile(datadir, 'savedir');
massdir = fullfile(datadir, 'masses'); 
marginaldir = fullfile(datadir, 'marginals');
lossdir = fullfile(datadir, 'pairwise_losses');
VOCdevkitpath = '/home/faruk/VOCdevkit'; %%% REPLACE HERE WITH PATH FOR VOCdevkit

if ~exist(divsoldir, 'dir')
    mkdir(divsoldir);
end
if ~exist(savedir, 'dir')
    mkdir(savedir);
end
if ~exist(massdir, 'dir')
    mkdir(massdir);
end
if ~exist(marginaldir, 'dir')
    mkdir(marginaldir);
end
if ~exist(lossdir, 'dir')
    mkdir(lossdir);
end

addpath(fullfile(VOCdevkitpath, 'VOCcode'));
addpath('./helper_functions');
addpath(genpath('./UGM_modified'));

fprintf('paths added...');

radius_range   = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50] %, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] %[0.4, 0.5, 0.7, 0.8, 0.9]  %[0.05, 0.1, 0.2, 0.5] % [0.01, 0.02] %[ 0, 0.01, 0.02] %[0.005, 0.01]
calT_range     = [0.1, 1, 10]; %[ 0.001, 0.01, 0.1, 1, 10] %[0.0001, 0.001, 0.008, 0.01, 0.05, 1, 5, 10, 50]
lambda_range   = [0.01, 0.05, 0.10]; %[ 0.01 0.05, 0.10, 0.15] %[0.01, 0.05, 0.10, 0.15] %[0.01, 0.05, 0.10, 0.15, 0.20]

NUM_WORKERS = 40;

if isempty(gcp('nocreate'))
     parpool(NUM_WORKERS);
end

overlap = 1;
noOverlap = 1 - overlap;
mstar = 0;

nim = 1449;
nummodes = 10;

dataset = 'VOC2012';
trainset = 'train';
testset =  'val';
id = 'comp6';
solcreator = 'o2p';

fprintf('Setting up initial stuff:\n');
load_prelim_data;

fprintf('Computing/loading pairwise losses...\n');
load_c3rf_data;

fprintf('Performing cross-validation...\n');
crossval_c3rf;
