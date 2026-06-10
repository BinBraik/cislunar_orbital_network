function [ix, iy, it] = atlas_bin_xyth(x, y, th, grid3)
%ATLAS_BIN_XYTH  Bin (x,y,theta) to voxel indices using grid edges.
% Out-of-bounds => NaN indices.
%
% Theta is wrapped to [-pi,pi) before binning.

thw = wrap_to_pi(th);

ix = discretize(x, grid3.x_edges);
iy = discretize(y, grid3.y_edges);
it = discretize(thw, grid3.th_edges);

end
