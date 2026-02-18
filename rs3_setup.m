function rs3_setup()
%RS3_SETUP  Add rs3 source folders to MATLAB path.
%
% Usage (recommended):
%   cd <repo-root>
%   rs3_setup
%
% After running this, you can execute runners in /scripts.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'scripts'));

fprintf('[rs3] Path configured (root + src + scripts).\n');
end
