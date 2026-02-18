
function rs3_cfg_validate(cfg)
%RS3_CFG_VALIDATE  Validate cfg struct schema and values (rs3).

req = {'sys','families','seed','grid','fan','propag','log','cache','refine','cand','par','io','diag'};
for k = 1:numel(req)
    if ~isfield(cfg, req{k})
        error('rs3:cfg:missingField', 'cfg.%s is missing', req{k});
    end
end

% --- System ---
mustPos(cfg.sys.RE_nd, 'cfg.sys.RE_nd');
mustPos(cfg.sys.RM_nd, 'cfg.sys.RM_nd');

% --- Seed / PO ---
mustPos(cfg.seed.Tf_scale, 'cfg.seed.Tf_scale');
mustIntGe(cfg.seed.N_dense, 2, 'cfg.seed.N_dense');
mustFinite(cfg.seed.y_eps, 'cfg.seed.y_eps');
mustPos(cfg.seed.ds_seed, 'cfg.seed.ds_seed');
mustIntGe(cfg.seed.minSegPts, 2, 'cfg.seed.minSegPts');

% --- Grid ---
mustPos(cfg.grid.Rdom, 'cfg.grid.Rdom');
mustPos(cfg.grid.dx, 'cfg.grid.dx');
mustPos(cfg.grid.dy, 'cfg.grid.dy');
mustPos(cfg.grid.dtheta, 'cfg.grid.dtheta');

% --- Fan ---
% Accept either dtheta_fan or dtheta (alias)
if isfield(cfg.fan,'dtheta_fan')
    mustPos(cfg.fan.dtheta_fan, 'cfg.fan.dtheta_fan');
elseif isfield(cfg.fan,'dtheta')
    mustPos(cfg.fan.dtheta, 'cfg.fan.dtheta');
else
    error('rs3:cfg:fan', 'cfg.fan must include dtheta_fan');
end
mustPos(cfg.fan.DV_cap_nd, 'cfg.fan.DV_cap_nd');

% --- Propagation ---
mustPos(cfg.propag.Tmax, 'cfg.propag.Tmax');
mustPos(cfg.propag.absTol, 'cfg.propag.absTol');
mustPos(cfg.propag.relTol, 'cfg.propag.relTol');
mustFinite(cfg.propag.v2tol, 'cfg.propag.v2tol');

% --- Logging ---
mustPos(cfg.log.step_len_factor, 'cfg.log.step_len_factor');
mustPos(cfg.log.maxstep_factor, 'cfg.log.maxstep_factor');
if ~isfield(cfg.log,'segwalk') || ~isstruct(cfg.log.segwalk)
    error('rs3:cfg:log', 'cfg.log.segwalk struct is required');
end
if ~islogical(cfg.log.segwalk.enable)
    error('rs3:cfg:log', 'cfg.log.segwalk.enable must be logical');
end
mustPos(cfg.log.segwalk.frac, 'cfg.log.segwalk.frac');

% --- Refinement ---
mustPos(cfg.refine.dx_min, 'cfg.refine.dx_min');
mustPos(cfg.refine.dy_min, 'cfg.refine.dy_min');
mustPos(cfg.refine.dtheta_min, 'cfg.refine.dtheta_min');
mustIntGe(cfg.refine.maxLevels, 0, 'cfg.refine.maxLevels');
mustIntGe(cfg.refine.maxRegions, 1, 'cfg.refine.maxRegions');
% --- Refinement 7 (new) – optional, no strict checks in Phase 0 ---
if isfield(cfg,'refine7')
    if ~isfield(cfg.refine7,'enable') || ~islogical(cfg.refine7.enable)
        warning('rs3:cfg:refine7', 'cfg.refine7.enable missing or not logical; using false.');
        cfg.refine7.enable = false;
    end
    % other fields are allowed to be missing (will use defaults in code)
end
% --- Candidate ---
mustIntGe(cfg.cand.K_per_voxel, 1, 'cfg.cand.K_per_voxel');
if isfield(cfg.cand, 'K_pairs_per_voxel')
    mustIntGe(cfg.cand.K_pairs_per_voxel, 1, 'cfg.cand.K_pairs_per_voxel');
end
mustIntGe(cfg.cand.maxPairsPerVoxel, 1, 'cfg.cand.maxPairsPerVoxel');

% --- Parallel ---
if ~islogical(cfg.par.enable)
    error('rs3:cfg:par', 'cfg.par.enable must be logical');
end
if ~strcmp(cfg.par.mode, 'jobs')
    error('rs3:cfg:par', 'cfg.par.mode must be ''jobs''.');
end

% --- IO ---
mustText(cfg.io.out_root, 'cfg.io.out_root');
mustText(cfg.io.fig_visible, 'cfg.io.fig_visible');

% --- Diagnostics ---
if ~islogical(cfg.diag.run_step3)
    error('rs3:cfg:diag', 'cfg.diag.run_step3 must be logical');
end
if ~islogical(cfg.diag.run_step4)
    error('rs3:cfg:diag', 'cfg.diag.run_step4 must be logical');
end

end

% ---------------- helpers ----------------
function mustPos(v, name)
if ~isscalar(v) || ~isnumeric(v) || ~isfinite(v) || v <= 0
    error('rs3:cfg:badValue', '%s must be a finite positive scalar', name);
end
end

function mustFinite(v, name)
if ~isscalar(v) || ~isnumeric(v) || ~isfinite(v)
    error('rs3:cfg:badValue', '%s must be a finite scalar', name);
end
end

function mustIntGe(v, lo, name)
if ~isscalar(v) || ~isnumeric(v) || ~isfinite(v) || floor(v)~=v || v < lo
    error('rs3:cfg:badValue', '%s must be an integer >= %d', name, lo);
end
end

function mustText(v, name)
if ~(ischar(v) || isstring(v))
    error('rs3:cfg:badValue', '%s must be a string', name);
end
end
