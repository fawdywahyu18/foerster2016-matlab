function [is_stable, max_eig, eig_vals] = check_mss(H, P)
% Check Mean Square Stability (MSS) criterion.
% Foerster et al. (2016), page 651.
%
% A solution is stable iff all eigenvalues of
%   (P' kron I_{nx^2}) * blkdiag(H_1 kron H_1, ..., H_ns kron H_ns)
% are inside the unit circle.

ns = length(H);
nx = size(H{1},1);
nx2 = nx^2;

blk = zeros(ns*nx2);
for i = 1:ns
    rows = (i-1)*nx2+1:i*nx2;
    blk(rows, rows) = kron(H{i}, H{i});
end

M = kron(P', eye(nx2)) * blk;
eig_vals = eig(M);
max_eig = max(abs(eig_vals));
is_stable = max_eig < 1;
end
