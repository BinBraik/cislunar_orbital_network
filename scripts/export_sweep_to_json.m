%% EXPORT_SWEEP_TO_JSON
% Generate a self-contained web/dv_map.html from sweep results.
%
% Reads:
%   atlas_sweep_results/sweep_DVmatrix_results.mat  — DV/TOF matrices (required)
%   atlas_network_results/network_results.mat        — centrality data (optional)
%
% If network_results.mat is absent the script computes harmonic closeness
% locally using the corrected formula from net_centrality.m Step 6.
%
% Writes:
%   web/dv_map.html   — fully self-contained, open directly in any browser

clear; clc;

% ── paths ─────────────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));

SWEEP_MAT   = fullfile(repoRoot, 'atlas_sweep_results', 'sweep_DVmatrix_results.mat');
NET_MAT     = fullfile(repoRoot, 'atlas_network_results', 'network_results.mat');
HTML_TMPL   = fullfile(repoRoot, 'web', 'dv_map_template.html');
HTML_OUT    = fullfile(repoRoot, 'web', 'dv_map.html');

if ~exist(fullfile(repoRoot, 'web'), 'dir')
    mkdir(fullfile(repoRoot, 'web'));
end

% ── load sweep data ────────────────────────────────────────────────────────────
fprintf('Loading sweep data: %s\n', SWEEP_MAT);
if ~isfile(SWEEP_MAT)
    error('export_sweep_to_json: sweep .mat not found:\n  %s\n', SWEEP_MAT);
end
S = load(SWEEP_MAT, 'DV_cap_list', 'Tmax_list', 'families', ...
         'DVmatrix_sweep', 'TOFmatrix_sweep');

DV_cap_list     = S.DV_cap_list(:);
Tmax_list       = S.Tmax_list(:);
families        = S.families(:);
DVmatrix_sweep  = S.DVmatrix_sweep;
TOFmatrix_sweep = S.TOFmatrix_sweep;

N     = numel(families);
nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);

shortNames = { ...
    'LyapL1'; 'LyapL2'; 'Cyc21'; 'Cyc11a'; 'Cyc11b'; 'Cyc32'; ...
    'R21S'; 'R21U'; 'R31S'; 'R31U'; 'R52S'; 'R52U'; 'DPO' };

if numel(shortNames) ~= N
    error('Short name list length (%d) != N (%d).', numel(shortNames), N);
end

% ── best case: element-wise min over all sweep cells ─────────────────────────
fprintf('Computing best-case matrix...\n');
bestDV  = inf(N, N);
bestTOF = nan(N, N);
for di = 1:nDV
    for dj = 1:nTmax
        M = DVmatrix_sweep{di, dj};
        T = TOFmatrix_sweep{di, dj};
        if isempty(M) || ~isnumeric(M), continue; end
        improved = isfinite(M) & (M < bestDV);
        bestDV(improved)  = M(improved);
        bestTOF(improved) = T(improved);
    end
end
bestDV(isinf(bestDV))   = NaN;
bestDV(1:N+1:end)       = NaN;   % diagonal → NaN
bestTOF(1:N+1:end)      = NaN;

% ── load or compute harmonic closeness ────────────────────────────────────────
hc_all    = NaN(N, nDV, nTmax);
is_ap_all = false(N, nDV, nTmax);

if isfile(NET_MAT)
    fprintf('Loading network centrality from: %s\n', NET_MAT);
    NET = load(NET_MAT);
    % Support both old field name (reach_all) and new (harmonic_closeness_all)
    if isfield(NET, 'harmonic_closeness_all') && ...
            isequal(size(NET.harmonic_closeness_all), [N, nDV, nTmax])
        hc_all = NET.harmonic_closeness_all;
    elseif isfield(NET, 'reach_all') && isequal(size(NET.reach_all), [N, nDV, nTmax])
        hc_all = NET.reach_all;
        warning('Using legacy reach_all field; re-run sweep for corrected harmonic_closeness_all.');
    else
        warning('harmonic_closeness_all not found or wrong size — recomputing locally.');
    end
    if isfield(NET, 'is_ap_all') && isequal(size(NET.is_ap_all), [N, nDV, nTmax])
        is_ap_all = NET.is_ap_all;
    end
else
    fprintf('network_results.mat not found — computing harmonic closeness locally...\n');
end

% Fill any remaining NaN snapshots by computing harmonic closeness from DV matrix
for di = 1:nDV
    for dj = 1:nTmax
        if all(isfinite(hc_all(:, di, dj))), continue; end
        M = DVmatrix_sweep{di, dj};
        if isempty(M) || ~isnumeric(M), continue; end
        hc_all(:, di, dj) = local_compute_harmonic_closeness(M, N);
    end
end

% best-case: use max-budget snapshot
bestReach = hc_all(:, nDV, nTmax);
bestIsAP  = is_ap_all(:, nDV, nTmax);

% ── build JSON ────────────────────────────────────────────────────────────────
fprintf('Building JSON...\n');

% Helper: convert NaN/Inf to null in JSON string (JSON spec has no NaN/Inf)
nanToNull = @(s) regexprep(s, '-?Inf(inity)?|-?NaN', 'null');

% Convert 2-D matrix to JSON 2-D array string
mat2json  = @(M) nanToNull(jsonencode(M));

% bestCase block
bestCaseStr = sprintf('{"dv":%s,"tof":%s,"reach":%s}', ...
    mat2json(bestDV), ...
    mat2json(bestTOF), ...
    nanToNull(jsonencode(bestReach(:)')));

% sweep array
sweepParts = cell(nDV * nTmax, 1);
k = 0;
for di = 1:nDV
    for dj = 1:nTmax
        k = k + 1;
        M  = DVmatrix_sweep{di, dj};
        T  = TOFmatrix_sweep{di, dj};
        R  = hc_all(:, di, dj);
        AP = is_ap_all(:, di, dj);

        if isempty(M) || ~isnumeric(M)
            dvStr  = 'null';
            tofStr = 'null';
        else
            M(1:N+1:end) = NaN;   % diagonal → NaN
            T(1:N+1:end) = NaN;
            dvStr  = mat2json(M);
            tofStr = mat2json(T);
        end
        sweepParts{k} = sprintf( ...
            '{"di":%d,"dj":%d,"dvCap":%.6g,"tmax":%.6g,"dv":%s,"tof":%s,"reach":%s,"isAP":%s}', ...
            di - 1, dj - 1, ...                  % 0-indexed for JS
            DV_cap_list(di), Tmax_list(dj), ...
            dvStr, tofStr, ...
            nanToNull(jsonencode(R(:)')), ...
            jsonencode(logical(AP(:)')));
    end
end

sweepStr = ['[', strjoin(sweepParts, ','), ']'];

% top-level DATA object
dataJson = sprintf([ ...
    '{\n' ...
    '  "families": %s,\n' ...
    '  "shortNames": %s,\n' ...
    '  "dvCapList": %s,\n' ...
    '  "tmaxList": %s,\n' ...
    '  "bestCase": %s,\n' ...
    '  "sweep": %s\n' ...
    '}'], ...
    jsonencode(families(:)'), ...
    jsonencode(shortNames(:)'), ...
    jsonencode(DV_cap_list(:)'), ...
    jsonencode(Tmax_list(:)'), ...
    bestCaseStr, ...
    sweepStr);

% ── read template and inject data ─────────────────────────────────────────────
if isfile(HTML_TMPL)
    fprintf('Reading template: %s\n', HTML_TMPL);
    fid = fopen(HTML_TMPL, 'r');
    htmlStr = fread(fid, '*char')';
    fclose(fid);
    htmlStr = strrep(htmlStr, '/*INJECT_DATA*/null/*END_INJECT*/', dataJson);
else
    % No template: embed data into the standalone HTML (built below)
    htmlStr = [];
end

% ── write output ──────────────────────────────────────────────────────────────
if isempty(htmlStr)
    % Template not found — just write data as JS for manual use
    outJs = fullfile(repoRoot, 'web', 'sweep_data.js');
    fid = fopen(outJs, 'w');
    fprintf(fid, 'const DATA = %s;\n', dataJson);
    fclose(fid);
    fprintf('\nNo template found. Data written to:\n  %s\n', outJs);
    fprintf('Place dv_map_template.html alongside it to generate dv_map.html.\n');
else
    fid = fopen(HTML_OUT, 'w');
    fwrite(fid, htmlStr);
    fclose(fid);
    fprintf('\nSelf-contained HTML written to:\n  %s\n', HTML_OUT);
    fprintf('Open it directly in Chrome or Firefox.\n');

    % Keep the published GitHub Pages copy (docs/index.html) in sync
    PAGES_OUT = fullfile(repoRoot, 'docs', 'index.html');
    if isfile(PAGES_OUT)
        copyfile(HTML_OUT, PAGES_OUT);
        fprintf('Also updated GitHub Pages copy:\n  %s\n', PAGES_OUT);
    end
end

fprintf('Done.\n');

% ── local functions ────────────────────────────────────────────────────────────

function reach = local_compute_harmonic_closeness(M, N)
%LOCAL_COMPUTE_HARMONIC_CLOSENESS  Harmonic closeness (net_centrality Step 6).
%   Uses Floyd-Warshall on the raw DV matrix; treats NaN as no direct edge.
%   Applies budget-reachability filter: pairs with dist > max(M(:)) contribute 0.

W = M;
W(isnan(W)) = Inf;
W(1:N+1:end) = 0;

% Floyd-Warshall
dist = W;
for kk = 1:N
    for ii = 1:N
        for jj = 1:N
            alt = dist(ii, kk) + dist(kk, jj);
            if alt < dist(ii, jj)
                dist(ii, jj) = alt;
            end
        end
    end
end

reach = zeros(N, 1);
for ii = 1:N
    finite_dests = find(isfinite(dist(ii, :)) & (1:N) ~= ii);
    if isempty(finite_dests), continue; end
    reach(ii) = sum(1 ./ dist(ii, finite_dests)) / (N - 1);
end
end
