function S = rs3_family_build_hits(S, cfg, grid3)
%RS3_FAMILY_BUILD_HITS  Step 4: propagate fan from seeds and log voxel hits.
%
% Builds only UPPER-half by integration. Lower-half is computed on-the-fly
% in Step 6 via rs3_rows_mirror_lower (Phase 4 memory optimization).
%
% CHANGES from original:
%   - Rows stored as packed struct-of-arrays (~4x less memory)
%   - Lower-half rows NOT stored (saves additional 50%)
%   - Aggressive clear after vertcat (Phase 2)
%   - rs3_rows_concat_cell replaces vertcat for packed structs
%
% Outputs added to S.Step4:
%   delta_lists, stats, rows_FRS_upper, rows_BRS_upper (packed structs)
%
% Row struct fields: iSeed, iHead, leg, t, ix, iy, it, halfFlag, n

tStart = tic;

SeedsU = S.SeedsUpper;
Nseeds = size(SeedsU,1);

% ---- Build per-seed delta lists (fan wedge) ----
[delta_lists, fanStats] = rs3_fan_delta_lists(SeedsU, S.CJ, S.mu, cfg);

% ---- Build job list (seed x heading) ----
jobSeed = zeros(0,1);
jobHead = zeros(0,1);
jobDelta = zeros(0,1);

for i = 1:Nseeds
    deltas = delta_lists{i};
    nh = numel(deltas);
    jobSeed = [jobSeed; i*ones(nh,1)]; %#ok<AGROW>
    jobHead = [jobHead; (1:nh)'];      %#ok<AGROW>
    jobDelta = [jobDelta; deltas(:)];  %#ok<AGROW>
end

nJobs = numel(jobSeed);
fprintf('[rs3] Step 4 (%s): jobs=%d (seeds=%d, mean heads=%.2f, max heads=%d)\n', ...
        S.name, nJobs, Nseeds, fanStats.nheads_mean, fanStats.nheads_max);

% ---- ODE options ----
RE = cfg.sys.RE_nd;
RM = cfg.sys.RM_nd;
step_len = cfg.log.step_len_factor * min(grid3.dx, grid3.dy);
MaxStep = cfg.log.maxstep_factor * step_len;

optsR = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol, 'MaxStep', MaxStep, ...
    'Events', @(t,X) rs3_ev_stop_reduced(t, X, S.mu, RE, RM, grid3.Rdom, S.CJ, cfg.propag.v2tol));
opts4 = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol, 'MaxStep', MaxStep, ...
    'Events', @(t,X) rs3_ev_stop_full4d(t, X, S.mu, RE, RM, grid3.Rdom));

Tmax = cfg.propag.Tmax;

% ---- Parfor over jobs (seed x heading) ----
usePar = cfg.par.enable;
if usePar && isempty(gcp('nocreate'))
    warning('[rs3] cfg.par.enable=true but no parpool exists. Running serial.');
    usePar = false;
end

rowsF_cell = cell(nJobs,1);
rowsB_cell = cell(nJobs,1);

if usePar
    fprintf('[rs3]   Step4: running PARFOR over seed×heading jobs.\n');

    doProg = cfg.diag.progress && isfield(cfg.par,'progress_every') && ~isempty(cfg.par.progress_every);
    dq = []; %#ok<NASGU>
    if doProg
        dq = parallel.pool.DataQueue;
        rs3_progress_tick('init', nJobs, cfg.par.progress_every);
        afterEach(dq, @(~) rs3_progress_tick('tick'));
    end

    parfor j = 1:nJobs
        iSeed = jobSeed(j);
        iHead = jobHead(j);
        delta = jobDelta(j);

        x0 = SeedsU(iSeed,1);
        y0 = SeedsU(iSeed,2);
        th0 = SeedsU(iSeed,3);
        th = rs3_wrapToPi(th0 + delta);

        X0 = [x0; y0; th];

        [tF, XF] = rs3_core_integrate_reduced_with_fallback(X0, [0 Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsF_cell{j} = rs3_hits_log_from_traj(iSeed, iHead, 1, +1, tF, XF, grid3, cfg);

        [tB, XB] = rs3_core_integrate_reduced_with_fallback(X0, [0 -Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsB_cell{j} = rs3_hits_log_from_traj(iSeed, iHead, 2, +1, tB, XB, grid3, cfg);
        if doProg
            send(dq, 1);
        end
    end
else
    tProg = tic;
    every = max(1, round(nJobs/20));
    for j = 1:nJobs
        iSeed = jobSeed(j);
        iHead = jobHead(j);
        delta = jobDelta(j);

        x0 = SeedsU(iSeed,1);
        y0 = SeedsU(iSeed,2);
        th0 = SeedsU(iSeed,3);
        th = rs3_wrapToPi(th0 + delta);

        X0 = [x0; y0; th];

        [tF, XF] = rs3_core_integrate_reduced_with_fallback(X0, [0 Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsF_cell{j} = rs3_hits_log_from_traj(iSeed, iHead, 1, +1, tF, XF, grid3, cfg);

        [tB, XB] = rs3_core_integrate_reduced_with_fallback(X0, [0 -Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsB_cell{j} = rs3_hits_log_from_traj(iSeed, iHead, 2, +1, tB, XB, grid3, cfg);

        if cfg.diag.progress && (mod(j, every)==0 || j==nJobs)
            fprintf('[rs3]   Step4 progress: %d/%d (%.1f%%) | elapsed %.1fs\n', ...
                j, nJobs, 100*j/nJobs, toc(tProg));
        end
    end
end

% ---- Concatenate packed structs (replaces vertcat) ----
rows_FRS_upper = rs3_rows_concat_cell(rowsF_cell);
clear rowsF_cell;                  % Phase 2: free cell array immediately
rows_BRS_upper = rs3_rows_concat_cell(rowsB_cell);
clear rowsB_cell;                  % Phase 2: free cell array immediately

% ---- Validate hits ----
rs3_hits_validate(rows_FRS_upper, rows_BRS_upper, grid3);

% ---- Phase 4: do NOT compute or store lower-half rows ----
% Lower-half is generated on-the-fly in Step 6 via rs3_rows_mirror_lower.
% This saves ~50% of permanent atlas memory.

S.Step4 = struct();
S.Step4.Tmax = Tmax;
S.Step4.step_len = step_len;
S.Step4.MaxStep = MaxStep;
S.Step4.delta_lists = delta_lists;
S.Step4.fanStats = fanStats;
S.Step4.rows_FRS_upper = rows_FRS_upper;
S.Step4.rows_BRS_upper = rows_BRS_upper;
S.Step4.nJobs = nJobs;
S.Step4.packed = true;  % flag for downstream compatibility checks

fprintf('[rs3] Step 4 (%s) done in %.2fs. rows(FRS_u)=%d rows(BRS_u)=%d\n', ...
        S.name, toc(tStart), double(rows_FRS_upper.n), double(rows_BRS_upper.n));
end
