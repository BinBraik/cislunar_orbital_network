function S = atlas_family_build_hits(S, cfg, grid3)
%RS3_FAMILY_BUILD_HITS  Step 4: propagate fan from seeds and log voxel hits.
%
% Exploits CR3BP time-reversal symmetry: BRS = R(FRS).
% Integrates FORWARD only, from both upper and lower seeds.
%   rows_FRS_upper : upper-half hits from forward trajectories of upper seeds
%   rows_FRS_lower : lower-half hits from forward trajectories of lower seeds
% BRS is derived on-the-fly downstream via R(FRS): no backward integration.
%
% CHANGES from original:
%   - Rows stored as packed struct-of-arrays (~4x less memory)
%   - Lower-half rows stored directly (forward from lower seeds)
%   - No backward integration: BRS = R(FRS) applied downstream
%   - Aggressive clear after vertcat (Phase 2)
%   - atlas_rows_concat_cell replaces vertcat for packed structs
%
% Outputs added to S.Step4:
%   delta_lists, stats, rows_FRS_upper, rows_FRS_lower (packed structs)
%
% Row struct fields: iSeed, iHead, leg, t, ix, iy, it, halfFlag, n

tStart = tic;

% Ensure per-family Keep mask (CJ + primary radii) is used during hit logging.
if isfield(S,'grid3') && isstruct(S.grid3) && isfield(S.grid3,'Keep')
    grid3 = S.grid3;
end

SeedsU = S.SeedsUpper;
SeedsL = S.SeedsLower;
Nseeds = size(SeedsU,1);

% ---- Build per-seed delta lists (fan wedge) ----
% Delta lists are identical for upper and lower seeds (same speed by symmetry).
[delta_lists, fanStats] = atlas_fan_delta_lists(SeedsU, S.CJ, S.mu, cfg);

% ---- Build job list: upper seeds then lower seeds ----
% iSeed is always in 1..Nseeds (same index used for both upper and lower).
% jobIsUpper distinguishes which seed matrix to look up.
jobSeed  = zeros(0,1);
jobHead  = zeros(0,1);
jobDelta = zeros(0,1);

for i = 1:Nseeds
    deltas = delta_lists{i};
    nh = numel(deltas);
    jobSeed  = [jobSeed;  i*ones(nh,1)]; %#ok<AGROW>
    jobHead  = [jobHead;  (1:nh)'];      %#ok<AGROW>
    jobDelta = [jobDelta; deltas(:)];    %#ok<AGROW>
end

nJobsHalf = numel(jobSeed);          % jobs per half (upper or lower)
nJobs     = 2 * nJobsHalf;           % total jobs (forward only, both halves)

% Concatenate: first nJobsHalf entries = upper seeds, next = lower seeds
jobSeed_all  = [jobSeed;  jobSeed];
jobHead_all  = [jobHead;  jobHead];
jobDelta_all = [jobDelta; jobDelta];
jobIsUpper   = [true(nJobsHalf,1); false(nJobsHalf,1)];

fprintf('[atlas] Step 4 (%s): jobs=%d (seeds=%d x2 halves, mean heads=%.2f, max heads=%d)\n', ...
        S.name, nJobs, Nseeds, fanStats.nheads_mean, fanStats.nheads_max);

% ---- ODE options ----
RE = cfg.sys.RE_nd;
RM = cfg.sys.RM_nd;
step_len = cfg.log.step_len_factor * min(grid3.dx, grid3.dy);
MaxStep = cfg.log.maxstep_factor * step_len;

optsR = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol, 'MaxStep', MaxStep, ...
    'Events', @(t,X) cr3bp_ev_stop_reduced(t, X, S.mu, RE, RM, grid3.Rdom, S.CJ, cfg.propag.v2tol));
opts4 = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol, 'MaxStep', MaxStep, ...
    'Events', @(t,X) cr3bp_ev_stop_full(t, X, S.mu, RE, RM, grid3.Rdom));

Tmax = cfg.propag.Tmax;

% ---- Parfor over jobs (seed x heading) ----
usePar = cfg.par.enable;
if usePar && isempty(gcp('nocreate'))
    warning('[atlas] cfg.par.enable=true but no parpool exists. Running serial.');
    usePar = false;
end

rowsF_cell = cell(nJobs,1);   % all forward hits (upper half: 1..nJobsHalf, lower: nJobsHalf+1..nJobs)

if usePar
    fprintf('[atlas]   Step4: running PARFOR over seed×heading jobs (forward only, both halves).\n');

    doProg = cfg.diag.progress && isfield(cfg.par,'progress_every') && ~isempty(cfg.par.progress_every);
    dq = []; %#ok<NASGU>
    if doProg
        dq = parallel.pool.DataQueue;
        progress_tick('init', nJobs, cfg.par.progress_every);
        afterEach(dq, @(~) progress_tick('tick'));
    end

    parfor j = 1:nJobs
        iSeed  = jobSeed_all(j);
        iHead  = jobHead_all(j);
        delta  = jobDelta_all(j);
        isUpper = jobIsUpper(j);

        if isUpper
            x0 = SeedsU(iSeed,1);
            y0 = SeedsU(iSeed,2);
            th0 = SeedsU(iSeed,3);
            hf  = int8(+1);
        else
            x0 = SeedsL(iSeed,1);
            y0 = SeedsL(iSeed,2);
            th0 = SeedsL(iSeed,3);
            hf  = int8(-1);
        end

        th = wrap_to_pi(th0 + delta);
        X0 = [x0; y0; th];

        [tF, XF] = cr3bp_integrate(X0, [0 Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsF_cell{j} = atlas_hits_log_from_traj(iSeed, iHead, 1, double(hf), tF, XF, grid3, cfg);

        if doProg
            send(dq, 1);
        end
    end
else
    tProg = tic;
    every = max(1, round(nJobs/20));
    for j = 1:nJobs
        iSeed  = jobSeed_all(j);
        iHead  = jobHead_all(j);
        delta  = jobDelta_all(j);
        isUpper = jobIsUpper(j);

        if isUpper
            x0 = SeedsU(iSeed,1);
            y0 = SeedsU(iSeed,2);
            th0 = SeedsU(iSeed,3);
            hf  = +1;
        else
            x0 = SeedsL(iSeed,1);
            y0 = SeedsL(iSeed,2);
            th0 = SeedsL(iSeed,3);
            hf  = -1;
        end

        th = wrap_to_pi(th0 + delta);
        X0 = [x0; y0; th];

        [tF, XF] = cr3bp_integrate(X0, [0 Tmax], S.CJ, S.mu, cfg.propag.v2tol, optsR, opts4);
        rowsF_cell{j} = atlas_hits_log_from_traj(iSeed, iHead, 1, hf, tF, XF, grid3, cfg);

        if cfg.diag.progress && (mod(j, every)==0 || j==nJobs)
            fprintf('[atlas]   Step4 progress: %d/%d (%.1f%%) | elapsed %.1fs\n', ...
                j, nJobs, 100*j/nJobs, toc(tProg));
        end
    end
end

% ---- Concatenate packed structs (replaces vertcat) ----
% First nJobsHalf cells = upper-seed forward hits (halfFlag=+1)
% Next  nJobsHalf cells = lower-seed forward hits (halfFlag=-1)
rows_FRS_upper = atlas_rows_concat_cell(rowsF_cell(1:nJobsHalf));
rows_FRS_lower = atlas_rows_concat_cell(rowsF_cell(nJobsHalf+1:end));
clear rowsF_cell;                  % Phase 2: free cell array immediately

% ---- Validate hits ----
atlas_hits_validate(rows_FRS_upper, rows_FRS_lower, grid3);

% BRS is NOT stored — derived on-the-fly as R(FRS) downstream.

S.Step4 = struct();
S.Step4.Tmax = Tmax;
S.Step4.step_len = step_len;
S.Step4.MaxStep = MaxStep;
S.Step4.delta_lists = delta_lists;
S.Step4.fanStats = fanStats;
S.Step4.rows_FRS_upper = rows_FRS_upper;
S.Step4.rows_FRS_lower = rows_FRS_lower;
S.Step4.nJobs = nJobs;
S.Step4.packed = true;  % flag for downstream compatibility checks

fprintf('[atlas] Step 4 (%s) done in %.2fs. rows(FRS_u)=%d rows(FRS_l)=%d\n', ...
        S.name, toc(tStart), double(rows_FRS_upper.n), double(rows_FRS_lower.n));
end
