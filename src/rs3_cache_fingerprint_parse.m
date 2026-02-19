function p = rs3_cache_fingerprint_parse(fp)
%RS3_CACHE_FINGERPRINT_PARSE  Parse fingerprint string into a readable struct.
%
% Input:
%   fp  — fingerprint string from cache_meta.fingerprint
%         Format (from rs3_cache_fingerprint_family):
%           rs3|{family}|mu=...|CJ=...|grid:R=...,dx=...,dy=...,dth=...|
%           seed:TfS=...,N=...,ds=...,yE=...,minSeg=...|fan:dth=...,DV=...|
%           prop:T=...,rt=...,at=...,v2=...|log:seg=...,frac=...,slf=...,msf=...|ver=...
%
% Output:
%   p   — struct with fields:
%           family, mu, CJ,
%           grid.{R, dx, dy, dtheta},
%           seed.{Tf_scale, N_dense, ds_seed, y_eps, minSegPts},
%           fan.{dtheta_fan, DV_cap_nd},
%           propag.{Tmax, relTol, absTol, v2tol},
%           log.{segwalk_enable, segwalk_frac, step_len_factor, maxstep_factor},
%           version_tag
%           raw   — original fingerprint string

p = struct();
p.raw = fp;
p.family = '';
p.mu = NaN; p.CJ = NaN;
p.grid   = struct('R',NaN,'dx',NaN,'dy',NaN,'dtheta',NaN);
p.seed   = struct('Tf_scale',NaN,'N_dense',NaN,'ds_seed',NaN,'y_eps',NaN,'minSegPts',NaN);
p.fan    = struct('dtheta_fan',NaN,'DV_cap_nd',NaN);
p.propag = struct('Tmax',NaN,'relTol',NaN,'absTol',NaN,'v2tol',NaN);
p.log    = struct('segwalk_enable',NaN,'segwalk_frac',NaN,'step_len_factor',NaN,'maxstep_factor',NaN);
p.version_tag = '';

if isempty(fp), return; end

parts = strsplit(fp, '|');
% parts{1} = 'rs3'
% parts{2} = familyName
% parts{3} = 'mu=...'
% parts{4} = 'CJ=...'
% parts{5} = 'grid:...'
% parts{6} = 'seed:...'
% parts{7} = 'fan:...'
% parts{8} = 'prop:...'
% parts{9} = 'log:...'
% parts{10} = 'ver=...'

if numel(parts) >= 2, p.family = parts{2}; end
if numel(parts) >= 3, p.mu = local_val(parts{3}, 'mu=');        end
if numel(parts) >= 4, p.CJ = local_val(parts{4}, 'CJ=');        end

if numel(parts) >= 5
    g = local_kv_block(parts{5}, 'grid:');
    p.grid.R      = local_kv(g, 'R');
    p.grid.dx     = local_kv(g, 'dx');
    p.grid.dy     = local_kv(g, 'dy');
    p.grid.dtheta = local_kv(g, 'dth');
end

if numel(parts) >= 6
    s = local_kv_block(parts{6}, 'seed:');
    p.seed.Tf_scale  = local_kv(s, 'TfS');
    p.seed.N_dense   = local_kv(s, 'N');
    p.seed.ds_seed   = local_kv(s, 'ds');
    p.seed.y_eps     = local_kv(s, 'yE');
    p.seed.minSegPts = local_kv(s, 'minSeg');
end

if numel(parts) >= 7
    f = local_kv_block(parts{7}, 'fan:');
    p.fan.dtheta_fan = local_kv(f, 'dth');
    p.fan.DV_cap_nd  = local_kv(f, 'DV');
end

if numel(parts) >= 8
    r = local_kv_block(parts{8}, 'prop:');
    p.propag.Tmax   = local_kv(r, 'T');
    p.propag.relTol = local_kv(r, 'rt');
    p.propag.absTol = local_kv(r, 'at');
    p.propag.v2tol  = local_kv(r, 'v2');
end

if numel(parts) >= 9
    l = local_kv_block(parts{9}, 'log:');
    p.log.segwalk_enable    = local_kv(l, 'seg');
    p.log.segwalk_frac      = local_kv(l, 'frac');
    p.log.step_len_factor   = local_kv(l, 'slf');
    p.log.maxstep_factor    = local_kv(l, 'msf');
end

if numel(parts) >= 10
    v = strtrim(parts{10});
    if strncmp(v, 'ver=', 4)
        p.version_tag = v(5:end);
    end
end
end

% ---- helpers ----
function v = local_val(str, prefix)
% Extract scalar value from 'prefix=value' string.
v = NaN;
if strncmp(str, prefix, numel(prefix))
    v = str2double(str(numel(prefix)+1:end));
end
end

function kv = local_kv_block(str, prefix)
% Strip prefix ('grid:' etc.) and return the remaining key=val,... string.
if strncmp(str, prefix, numel(prefix))
    kv = str(numel(prefix)+1:end);
else
    kv = str;
end
end

function v = local_kv(kvstr, key)
% Extract value of 'key=value' from a comma-separated kv string.
v = NaN;
toks = strsplit(kvstr, ',');
prefix = [key '='];
for i = 1:numel(toks)
    t = strtrim(toks{i});
    if strncmp(t, prefix, numel(prefix))
        v = str2double(t(numel(prefix)+1:end));
        return;
    end
end
end
