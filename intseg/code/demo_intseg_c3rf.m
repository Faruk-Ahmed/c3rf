%% add paths:
clear; close; clc;
addpath(genpath('./'));

%% initialize paths:
datadir = '../data/voctest50data';
divsoldir = '../data/divsols';
gtdir = fullfile(datadir, 'gtdir');
savedir_masses = '../data/masses';
savedir_marginals = '../data/marginals';
savedir_pairwise_losses = '../data/pairwise_losses';

if ~exist(savedir_masses, 'dir')
    mkdir(savedir_masses);
end
if ~exist(savedir_marginals, 'dir')
    mkdir(savedir_marginals);
end
if ~exist(savedir_pairwise_losses, 'dir')
    mkdir(savedir_pairwise_losses);
end

flist = dir(fullfile(datadir,'*.mat'));
numimages = 50;
nummodes = 30;

%% Set parameter grids:
lambda_range = [0.1, 0.2, 0.3];  %%% [0.1, 0.2, 0.3, 0.4];
calT_range   = [1, 10, 100, 1000];    %%% [0.01, 0.1, 1, 10, 100, 1000, 100000];
radius_range = [0, 0.01, 0.05, 0.1]; %%% [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
noOverlap    = 0;
%% load/compute image data and masses and marginals:
load_prelim_data;

%% load/compute expected losses:
load_c3rf_data;

%% perform prediction:
crossval_c3rf;
