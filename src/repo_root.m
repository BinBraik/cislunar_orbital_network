function root = repo_root()
%REPO_ROOT  Return the repository root directory (folder that contains /src).
%
% This helper keeps paths stable regardless of the current working directory.

thisFile = mfilename('fullpath');
root = fileparts(fileparts(thisFile));
end
