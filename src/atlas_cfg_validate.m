
function atlas_cfg_validate(cfg)
%RS3_CFG_VALIDATE  Validate cfg struct schema and values (rs3).

req = {'sys','families','seed','grid','fan','propag','log','cache','par','io','diag'};
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
mustPos(cfg.grid.Rdom,   'cfg.grid.Rdom');
mustPos(cfg.grid.dx,     'cfg.grid.dx');
mustPos(cfg.grid.dy,     'cfg.grid.dy');
mustPos(cfg.grid.dtheta, 'cfg.grid.dtheta');

% --- Fan ---
mustPos(cfg.fan.dtheta_fan, 'cfg.fan.dtheta_fan');
mustPos(cfg.fan.DV_cap_nd,  'cfg.fan.DV_cap_nd');

% --- Propagation ---
mustPos(cfg.propag.Tmax,   'cfg.propag.Tmax');
mustPos(cfg.propag.absTol, 'cfg.propag.absTol');
mustPos(cfg.propag.relTol, 'cfg.propag.relTol');
mustFinite(cfg.propag.v2tol, 'cfg.propag.v2tol');

% --- Logging ---
mustPos(cfg.log.step_len_factor, 'cfg.log.step_len_factor');
mustPos(cfg.log.maxstep_factor,  'cfg.log.maxstep_factor');
if ~isfield(cfg.log,'segwalk') || ~isstruct(cfg.log.segwalk)
    error('rs3:cfg:log', 'cfg.log.segwalk struct is required');
end
if ~islogical(cfg.log.segwalk.enable)
    error('rs3:cfg:log', 'cfg.log.segwalk.enable must be logical');
end
mustPos(cfg.log.segwalk.frac, 'cfg.log.segwalk.frac');

% --- Parallel ---
if ~islogical(cfg.par.enable)
    error('rs3:cfg:par', 'cfg.par.enable must be logical');
end

% --- IO ---
mustText(cfg.io.out_root,    'cfg.io.out_root');
mustText(cfg.io.fig_visible, 'cfg.io.fig_visible');

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
