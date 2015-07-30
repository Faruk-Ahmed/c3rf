clear; clc; close;
% The MATLAB package 'boundedline' should be downloaded and be in the path. Get it here:
% http://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

% paste matrices of num_folds X num_divsols values for the variables below (or load from mat file, whatever)
D = % Delta-EIoEU+enum      (Run Premachandran et al.)
E = % C3RF-EIoEU+enum
M = % Mass-EIoEU+enum       (Run Premachandran et al. substituting masses in place of scores.)
F = % CRF-EIoEU+enum
MAP = % MAP
N = % CRF-EIoEU+greedy. You'll have to download Sebastian's code from 
    % http://research.microsoft.com/en-us/downloads/e164fe21-ef2b-4e34-98c1-4868968abb06/
    % and run it yourself.

NUM_SOLS  = 10;
NUM_FOLDS = 5; % this is actually, the number of permutations.

%%
clc; close;
[h, p] = boundedline(1:NUM_SOLS, mean(D), sqrt(var(D)/NUM_FOLDS), '-b*', ...
                     1:NUM_SOLS, mean(M), sqrt(var(M)/NUM_FOLDS), '-gO',  ...
                     1:NUM_SOLS, mean(E), sqrt(var(E)/NUM_FOLDS), '-m',  ...
                     1:NUM_SOLS, mean(F), sqrt(var(F)/NUM_FOLDS), '-c', ...
                     1:NUM_SOLS, mean(N), sqrt(var(N)/NUM_FOLDS), '--r', ...
                     1:NUM_SOLS, mean(MAP), sqrt(var(MAP)/NUM_FOLDS), '--k', 'alpha');

fsize = 14;
set(h, 'LineWidth', 2)
set(gcf, 'Color', 'w');
set(gca, 'YGrid', 'on')
set(gca, 'FontSize', fsize);

xlabel('#Solutions')
ylabel('Intersection-Over-Union')
title('Binary Foreground-Background Segmentation')

h1 = legend(h, {'Delta-EIoU+enum', 'Mass-EIoU+enum', 'C^3RF-EIoEU+enum', 'CRF-EIoEU+enum', 'CRF-EIoEU+greedy', 'MAP'}, 'Location', 'NorthWest');
set(h1 , 'FontSize', 12);
set(h1, 'interpreter', 'tex');
set(h1, 'Box', 'off');


