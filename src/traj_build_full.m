function traj = traj_build_full(T)
%RS5_BUILD_FULL_TRAJ  Build a single concatenated trajectory struct.
%
% Combines the departure arc A (forward from PO_A to patch) and the BRS
% arrival arc B (reversed, so it runs chronologically from patch toward
% PO_B) into a single trajectory that can be plotted or saved as a unit.
%
% Works with both:
%   T  — uncorrected struct from overlap_voxel_traj_extract
%   Tc — corrected  struct from traj_diffcorr
%
% The arc B stored in T/Tc is the R-transform of an FRS arc from PO_B:
%   x_B(1), y_B(1)  ≈ start near R-transformed PO_B
%   x_B(end), y_B(end) ≈ patch point
% Reversing it gives the physical path from patch toward PO_B_R.
%
% Output fields:
%   traj.x, .y, .th   — concatenated position/heading  [n_A+n_B × 1]
%   traj.t            — time   (arc A: 0..t_A;  arc B reversed: t_A..t_A+t_B)
%   traj.n_A          — number of points in arc A segment
%   traj.n_B          — number of points in arc B segment
%   traj.depart_x/y   — departure point (first point of arc A, on PO_A)
%   traj.depart_th    — heading at departure (IC_A(3) = seed_th + delta_A)
%   traj.patch_x/y    — patch / transfer point
%   traj.arrive_x/y   — arrival point  (last point of reversed arc B, near PO_B_R)
%   traj.arrive_th    — velocity direction at arrival (th_B(1) + pi)
%   traj.DV_*         — delta-V budget fields  (m/s)
%   traj.tof_*        — time-of-flight fields  (days)

% ---- extract common fields (present in both T and Tc) ----
XA     = T.XA;
tA_vec = T.tA_vec;
x_B    = T.x_B;
y_B    = T.y_B;
th_B   = T.th_B;
tB_vec = T.tB_vec;
IC_A   = T.IC_A;
t_A    = T.t_A;
t_B    = T.t_B;

% ---- patch point and total DV (field names differ between T and Tc) ----
if isfield(T, 'xp')
    % Corrected Tc: exact patch point from constraint satisfaction
    xpatch   = T.xp;
    ypatch   = T.yp;
    DV_total = T.DV_total_mps;
else
    % Uncorrected T: use arc A endpoint as approximate patch
    xpatch   = XA(end, 1);
    ypatch   = XA(end, 2);
    DV_total = T.DV_total_true_mps;
end

% ---- reverse arc B: chronologically from patch toward target orbit ----
idx_rev  = numel(x_B) : -1 : 1;
x_B_rev  = x_B(idx_rev);
y_B_rev  = y_B(idx_rev);
th_B_rev = th_B(idx_rev);
% Shift reversed time so arc B starts exactly at t_A
t_B_rev  = t_A + (t_B - tB_vec(idx_rev));

% ---- concatenate ----
traj.x  = [XA(:,1);  x_B_rev];
traj.y  = [XA(:,2);  y_B_rev];
traj.th = [XA(:,3);  th_B_rev];
traj.t  = [tA_vec;   t_B_rev];

traj.n_A = size(XA, 1);
traj.n_B = numel(x_B_rev);

% ---- key geometric points ----
traj.depart_x  = XA(1, 1);
traj.depart_y  = XA(1, 2);
traj.depart_th = IC_A(3);        % post-kick heading at departure

traj.patch_x   = xpatch;
traj.patch_y   = ypatch;

traj.arrive_x  = x_B_rev(end);   % last point of reversed B = first of stored B
traj.arrive_y  = y_B_rev(end);
% Physical arrival direction: the spacecraft approaches PO_B with velocity
% opposite to the BRS arc's initial heading (th_B(1) is the heading at the
% start of the stored arc in BRS frame; reversing time flips the direction).
traj.arrive_th = wrap_to_pi(th_B(1) + pi);

% ---- DV budget ----
traj.DV_turn_A_mps = T.DV_turn_A_mps;
traj.DV_patch_mps  = T.DV_patch_mps;
traj.DV_turn_B_mps = T.DV_turn_B_mps;
traj.DV_total_mps  = DV_total;

% ---- time of flight ----
traj.tof_A_days     = T.tof_A_days;
traj.tof_B_days     = T.tof_B_days;
traj.tof_total_days = T.tof_A_days + T.tof_B_days;
end
