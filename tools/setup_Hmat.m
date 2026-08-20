% ---- Software configuration
% Make hmmvp available
addpath 'hmmvp0.16';

% Declare global variable
global gamb;
gamb.use = 1;  % General flag for using certain features (context-specific)

% Optional H-matrix compression of the BEM matrix
gamb.hmat.use = 1;  % Use H-matrix compression as implemented in hmmvp
gamb.hmat.rerr = 1e-3;  % Element-wise relative error of approximation (~RelTol in ODE integration)

% Base directory for H-matrix files
dr = './hmat';  % Ensure this is relative to the current working directory


% Add the appropriate subdirectory to the MATLAB path based on block_id
if any(strcmpi(block_id, {'visco_block_I', 'Block_I'})) % Panda & Lindsey, 2024
    hsub = 'hmat_I';
elseif any(strcmpi(block_id, {'visco_block_II', 'Block_II'})) % Panda & Lindsey, 2024
    hsub = 'hmat_II';
elseif any(strcmpi(block_id, {'visco_block_III', 'Block_III'})) % Panda & Lindsey, 2024
    hsub = 'hmat_III';
elseif any(strcmpi(block_id, {'himalaya_visco', 'lindsey_2018'})) % Lindsey et al., 2018
    hsub = 'hmat_lindsey';
else
    error('Unknown block name: %s', block_id);
end

if ~exist(dr, 'dir')
    mkdir(dr);
end
if ~exist(fullfile(dr, hsub), 'dir')
    mkdir(fullfile(dr, hsub));
end

addpath(fullfile(dr, hsub));  % Matches with actual directory name

% % Define ne (number of elements) - assumes el is already defined
% % Replace this with the actual data loading if el isn't defined yet
% % Example: load('your_data_file.mat'); % Uncomment and adjust if needed
% if ~exist('el', 'var')
%     error('Variable el is not defined.');
% end

ne = size(el, 1);  % Ensure el is a valid matrix/array from your data

% Define fn_dec for the H-matrix filename
fn_dec = sprintf('ne%dhmat%drerr%3.2f', ne, gamb.hmat.use, log10(gamb.hmat.rerr));

% Set the H-matrix file path based on block_id
gamb.hmat.savefn = sprintf('%s/%s/Hmat_%s.dat', dr, hsub, fn_dec);

% Debugging: Display the H-matrix file path
disp(['H-matrix file path set to: ', gamb.hmat.savefn]);

% Compute new H-matrix if requested (assumes compute_hmat is defined)
% Set compute_hmat = true to generate the file, false to reuse it
if ~exist('compute_hmat', 'var')
    compute_hmat = false;
end
if compute_hmat
    CalcG_rake('Make', el, nd, rakes, gamb.hmat.rerr, gamb.hmat.savefn);
end

% Verify the file exists before proceeding
if ~isfile(gamb.hmat.savefn)
    error('H-matrix file not found at: %s. Set compute_hmat = true to generate it.', gamb.hmat.savefn);
end

% Initialize the H-matrix
[gamb.hmat.id, nnz_hmat] = hm_mvp('init', gamb.hmat.savefn, 4);

% ---- End of Software code
