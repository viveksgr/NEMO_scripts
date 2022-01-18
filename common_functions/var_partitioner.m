function [a_,b_,ab,t] = var_partitioner(a,b,t,fig_handle)
% Takes mutually correlated column vectors a,b and yields a_,b_,ab such that
% a = a_ +ab and b = b_+ab and a_,b_ and ab are mutually orthogonal
% Note that the scales for a_,b_ and ab will change, renormalize them if needed.
%--------------------------------------------------------------------------
% Vivek Sagar VivekSagar2016@u.northwestern.edu November 10, 2020
%--------------------------------------------------------------------------

a = a-mean(a);
a = a/norm(a);
b = b-mean(b);
b = b/norm(b);

if nargin<4
    fig_handle = false;
end

if nargin<3
% Normal vector choice
t = rand(size(a));
t = t-mean(t);
t = t/norm(t);
end

if (sum(isnan(a))+sum(isnan(b))>0)
    error('Vectors contain NaNs')
end

% Normalize the vectors and save normalization parameters
a_mean = mean(a);
a_norm = norm(a-a_mean);
b_mean = mean(b);
b_norm = norm(b-b_mean);
a = (a-a_mean)./a_norm;
b = (b-b_mean)./b_norm;

% Magnitude of mutual vector
ab_mod = sqrt((a'*b));

% All single letter symbols other than d are for vectors 
% Direction of mutual vector x: m.x = d
m = a+b;
d = 2*ab_mod;

% Choose some normal n to m;
m_n = m/norm(m);
n = m_n*(m_n'*t)-t; % a is some vector not orthogonal to m
n = n/norm(n);

% x = d/|m| m_hat + sqrt(1-d^2/m^2) n
x = (d/norm(m))*m_n+sqrt(1-(d^2/(norm(m)^2)))*n;
ab = (ab_mod/norm(x))*x;

a_ = a-ab;
b_ = b-ab;
% a_ = a_*a_norm+a_mean;
% b_ = b_*b_norm+b_mean;

if fig_handle
    figure()
    corrs_ = [];
    corrs_(1) = fastcorr(a,b);
    corrs_(2) = fastcorr(a_,b_);
    corrs_(3) = fastcorr(a,a_);
    corrs_(4) = fastcorr(b,b_);
    corrs_(5) = fastcorr(a_,b);
    corrs_(6) = fastcorr(b_,a);   
    corrs_(7) = fastcorr(a_+b_,ab);
    bar(categorical({'1 a b','2 a_ b_','3 a a_','4 b b_','5 b a_', '6 a b_','7 ab , sum a_ b_'}),corrs_);
    savefig(fullfile(pwd,'varpartition.fig'))
end