function [delta_lists, stats] = rs3_fan_delta_lists(SeedsUpper, CJ, mu, cfg)
%RS3_FAN_DELTA_LISTS  Build per-seed heading delta lists (fan wedge).
%
% For each seed, compute speed v0 from CJ and U(x,y), then compute max heading
% change dmax allowed by cfg.fan.DV_cap_nd:
%   dv_turn = 2*v0*sin(|delta|/2) <= DVcap
%   => dmax = 2*asin( min(1, DVcap/(2*v0)) )
%
% Outputs
%   delta_lists : cell{Nseeds} with column vectors of deltas (rad), including endpoints
%   stats       : struct with counts, dmax stats

assert(isfield(cfg,'fan') && isstruct(cfg.fan), 'cfg.fan missing');
DVcap = cfg.fan.DV_cap_nd;
dtheta_fan = cfg.fan.dtheta_fan;

N = size(SeedsUpper,1);
delta_lists = cell(N,1);

dmax_all = zeros(N,1);
nheads = zeros(N,1);

for i = 1:N
    x0 = SeedsUpper(i,1); y0 = SeedsUpper(i,2);
    pot0 = rs3_core_cr3bp_U_and_derivs(x0, y0, mu);
    v0 = sqrt(max(2*pot0.U - CJ, 0));
    dmax = 2*asin(min(1, DVcap/(2*max(v0,eps))));
    dmax_all(i) = dmax;

    % Symmetric fan centered on delta=0 (nominal heading).
    % Nh is always odd so delta=0 is exactly the middle element,
    % and both endpoints ±dmax are exact (no floating-point drift).
    Nh = 2*floor(dmax/dtheta_fan) + 1;
    deltas = linspace(-dmax, dmax, Nh);
    delta_lists{i} = deltas(:);
    nheads(i) = numel(delta_lists{i});
end

stats = struct();
stats.Nseeds = N;
stats.nheads = nheads;
stats.nheads_mean = mean(nheads);
stats.nheads_max = max(nheads);
stats.dmax_min = min(dmax_all);
stats.dmax_max = max(dmax_all);
stats.dmax_mean = mean(dmax_all);
end
