function [num_r, num_c] = break_it(n)
ceil_ = ceil(sqrt(n));
floor_ = floor(sqrt(n));
num_r = ceil_;
if n<=ceil_*floor_
    num_c = floor_;
else
    num_c = ceil_;
end