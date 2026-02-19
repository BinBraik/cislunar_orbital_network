function S = rs3_family_prepare_seeds(familyName, cfg, grid3geom)
%RS3_FAMILY_PREPARE_SEEDS  Step 3: integrate PO, select seeds, build per-family Keep mask.
%
% Inputs
%   familyName : char/string
%   cfg        : rs3 config
%   grid3geom  : output of rs3_grid_make(cfg) (geometry only)
%
% Output
%   S : struct with:
%       .name, .mu, .CJ, .Tf_base, .Tf_PO, .X0, .Xpo
%       .SeedsUpper, .SeedsLower, .seedMeta
%       .grid3 (geometry + Keep)

% --- load family initial condition (allowed core reuse) ---
[mu, CJ, Tf_base, X0] = rs3_core_family_ic(char(familyName));
Tf_PO = cfg.seed.Tf_scale * Tf_base;

if cfg.io.verbose
    fprintf('    mu=%.15f | CJ=%.10f | Tf_base=%.6f | Tf_PO=%.6f\n', mu, CJ, Tf_base, Tf_PO);
end

% --- dense periodic orbit integration in reduced model [x;y;theta] ---
opts = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol);
solPO = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t, X, CJ, mu, true), [0 Tf_PO], X0, opts);

t_dense = linspace(0, Tf_PO, cfg.seed.N_dense);
Xpo = deval(solPO, t_dense)';

% --- seeds on ALL upper segments (y >= y_eps) ---
[SeedsUpper, SeedsLower, seedMeta] = rs3_seed_upper_segments_by_arclength(Xpo(:,1), ...
    Xpo(:,2), Xpo(:,3), cfg.seed.y_eps, cfg.seed.ds_seed, cfg.seed.minSegPts);

% --- per-family Keep mask (CJ dependent) ---
grid3 = grid3geom;
grid3.Keep = rs3_keep_mask_xy(grid3, CJ, mu, cfg.sys.RE_nd, cfg.sys.RM_nd);

% --- validations / sanity checks ---
if cfg.diag.smoke_checks
    % 1) upper-half enforcement
    if ~isempty(SeedsUpper)
        ymin = min(SeedsUpper(:,2));
        if ymin < cfg.seed.y_eps - 1e-12
            error('Step3: SeedsUpper violate y_eps: min(y)=%.3g < y_eps=%.3g', ymin, cfg.seed.y_eps);
        end
    end

    % 2) symmetry mapping checks (lower should match mirror)
    if ~isempty(SeedsUpper)
        yerr = max(abs(SeedsLower(:,2) + SeedsUpper(:,2)));
        th_mirror = rs3_wrapToPi(pi - SeedsUpper(:,3));
        therr = max(abs(rs3_circ_diff(SeedsLower(:,3), th_mirror)));
        if yerr > 1e-12 || therr > 1e-10
            error('Step3: symmetry mismatch: yerr=%.3g, therr=%.3g', yerr, therr);
        end
    end

    % 3) Keep mask should contain some area
    keepFrac = nnz(grid3.Keep) / numel(grid3.Keep);
    if keepFrac < 0.01
        warning('Step3: Keep mask is very small (%.3g). Check CJ/Rdom.', keepFrac);
    end
end

% --- pack ---
S = struct();
S.name = char(familyName);
S.mu = mu;
S.CJ = CJ;
S.Tf_base = Tf_base;
S.Tf_PO = Tf_PO;
S.X0 = X0(:);
S.t_dense = t_dense(:);
S.Xpo = Xpo;
S.SeedsUpper = SeedsUpper;
S.SeedsLower = SeedsLower;
S.seedMeta = seedMeta;
S.grid3 = grid3;

end
