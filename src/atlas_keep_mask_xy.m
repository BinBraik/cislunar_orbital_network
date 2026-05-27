function Keep = atlas_keep_mask_xy(grid3, CJ, mu, RE_nd, RM_nd)
%RS3_KEEP_MASK_XY  Allowed-mask on (x,y) cell centers.
% Keep(iy,ix)=true when ALL are satisfied at the cell center:
%   1) inside analysis domain (r <= Rdom)
%   2) energetically allowed (v^2 = 2U - CJ >= 0)
%   3) outside Earth radius (rE > RE_nd)
%   4) outside Moon radius  (rM > RM_nd)
%
% Notes:
% - Per-family, because CJ and mu differ across families.
% - RE_nd/RM_nd are optional for backward compatibility.

if nargin < 4 || isempty(RE_nd), RE_nd = 0; end
if nargin < 5 || isempty(RM_nd), RM_nd = 0; end

[Xc,Yc] = meshgrid(grid3.x_centers, grid3.y_centers);

inR = hypot(Xc,Yc) <= grid3.Rdom;

x1 = Xc + mu;      y1 = Yc;          % Earth-centered frame shift
x2 = Xc - 1 + mu;  y2 = Yc;          % Moon-centered frame shift

r1 = sqrt(x1.^2 + y1.^2);            % distance to Earth center
r2 = sqrt(x2.^2 + y2.^2);            % distance to Moon center

U = 0.5*(Xc.^2 + Yc.^2) + (1-mu)./r1 + mu./r2;
v2 = 2*U - CJ;

outsideEarth = r1 > RE_nd;
outsideMoon  = r2 > RM_nd;

Keep = inR & (v2 >= 0) & outsideEarth & outsideMoon;
end
