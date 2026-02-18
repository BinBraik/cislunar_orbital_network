function rs3_progress_tick(cmd, total, every, tag)
%RS3_PROGRESS_TICK  Lightweight progress ticker for PARFOR jobs.
%
% Usage:
%   rs3_progress_tick('init', nJobs, everyN)
%   rs3_progress_tick('init', nJobs, everyN, tagString)
%   rs3_progress_tick('tick')
%
% Designed to be used with parallel.pool.DataQueue + afterEach.
%
% Notes:
% - tagString is printed after the "[rs3]" prefix, useful to distinguish
%   Step4 vs Step7 refinement progress.

persistent done nTot nEvery t0 lastPrint label

if nargin < 1 || isempty(cmd)
    cmd = 'tick';
end

switch lower(cmd)
    case 'init'
        done = 0;
        nTot = total;
        nEvery = max(1, every);
        t0 = tic;
        lastPrint = 0;

        if nargin >= 4 && ~isempty(tag)
            label = tag;
        else
            label = 'Step4 progress';
        end

    case 'tick'
        done = done + 1;

        if done == nTot || (done - lastPrint) >= nEvery
            fprintf('[rs3]   %s: %d/%d (%.1f%%) | elapsed %.1fs\n', ...
                label, done, nTot, 100*done/nTot, toc(t0));
            lastPrint = done;
        end

    otherwise
        error('rs3_progress_tick:UnknownCmd', 'Unknown cmd: %s', cmd);
end

end
