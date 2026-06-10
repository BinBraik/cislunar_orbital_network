function setup()
%SETUP  Add source folders to MATLAB path.
%
% Usage (recommended):
%   cd <repo-root>
%   setup
%
% After running this, you can execute any runner in /scripts.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'src', 'network'));
addpath(fullfile(root, 'scripts'));

fprintf('[atlas] Path configured (root + src + src/network + scripts).\n');
end
