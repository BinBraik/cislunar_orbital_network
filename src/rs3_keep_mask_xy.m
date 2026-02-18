function Keep = rs3_keep_mask_xy(grid3, CJ, mu)
%RS3_KEEP_MASK_XY  CJ-allowed mask on (x,y) cell centers.
% Keep(iy,ix)=true if cell center is inside domain and v^2 = 2U - CJ >= 0.
%
% Notes:
% - Per-family, because CJ differs slightly across families.
% - Same intent as baseline compute_keep_mask.

[Xc,Yc] = meshgrid(grid3.x_centers, grid3.y_centers);

inR = hypot(Xc,Yc) <= grid3.Rdom;

x1 = Xc + mu;  y1 = Yc;
x2 = Xc - 1 + mu;  y2 = Yc;

r1 = sqrt(x1.^2 + y1.^2);
r2 = sqrt(x2.^2 + y2.^2);

U = 0.5*(Xc.^2 + Yc.^2) + (1-mu)./r1 + mu./r2;
v2 = 2*U - CJ;

Keep = inR & (v2 >= 0);
end
