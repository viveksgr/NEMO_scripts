function [max_curv, curvature] = find_curvature(curve)
% Given Nx1 curve, find the curvature function and the maximum curvature
% Assume uniform spacing. Use interp1 otherwise, first.
x = 1:1:length(curve);
y = curve;

dx  = gradient(x);
ddx = gradient(dx);
dy  = gradient(y);
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = (sqrt(dx .^2 + dy .^2)).^3;
curvature = abs(num ./ denom);
max_curv = max(abs(curvature));