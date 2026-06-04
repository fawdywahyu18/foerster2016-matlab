function Q = solve_rbc_second_order(sol, par, P, mu_vec, rho_vec, sigma_vec, method)
% Solve for second-order coefficients via numerical Hessian approach.
%
% The second-order policy function is:
%   y_hat(z) = A^y_i * z + 0.5 * z' * Q{i,1} * z
%   x_hat_k(z) = A^x_i(k,:) * z + 0.5 * z' * Q{i,k+1} * z
%
% where z = [xhat_{t-1}; eps_t; chi].
%
% At second order, D^2 F_i(0) / (dz_{l1} dz_{l2}) = 0 for all l1, l2.
% This is a linear system in the Q elements.

ns = size(P,1); ny = 1; nx = 2; n = ny + nx; nz = nx + 1 + 1;

gh_nodes   = [-sqrt(3), 0, sqrt(3)];
gh_weights = [1/6, 2/3, 1/6];

n_unique = nz*(nz+1)/2;
n_unknowns = ns * n * n_unique;

Q0 = cell(ns, n);
for i = 1:ns
    for v = 1:n
        Q0{i,v} = zeros(nz);
    end
end

% Compute forcing: Hessian of F with Q = 0
K = compute_hessian_vector(Q0, sol, par, P, mu_vec, rho_vec, sigma_vec, method, gh_nodes, gh_weights);

% Build linear operator column by column
L = zeros(n_unknowns, n_unknowns);
delta_val = 1.0;
for m = 1:n_unknowns
    Qp = Q0;
    [ri, vi, l1, l2] = linear_idx_to_Q(m, ns, n, nz);
    Qp{ri,vi}(l1,l2) = delta_val;
    Qp{ri,vi}(l2,l1) = delta_val;

    Km = compute_hessian_vector(Qp, sol, par, P, mu_vec, rho_vec, sigma_vec, method, gh_nodes, gh_weights);
    L(:,m) = (Km - K) / delta_val;
end

q_vec = L \ (-K);

Q = Q0;
for m = 1:n_unknowns
    [ri, vi, l1, l2] = linear_idx_to_Q(m, ns, n, nz);
    Q{ri,vi}(l1,l2) = q_vec(m);
    Q{ri,vi}(l2,l1) = q_vec(m);
end
end


function hvec = compute_hessian_vector(Q, sol, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w)
% Compute the vectorized Hessian of F_i(z) at z=0 for all regimes.
% Returns a vector of length ns * n * n_unique.

ns = size(P,1); n = 3; nz = 4;
n_unique = nz*(nz+1)/2;
h = 1e-4;

F0 = zeros(n, ns);
for i = 1:ns
    F0(:,i) = eval_rbc_F(zeros(nz,1), i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
end

hvec = zeros(ns * n * n_unique, 1);

for i = 1:ns
    m_idx = 0;
    for l1 = 1:nz
        for l2 = l1:nz
            m_idx = m_idx + 1;

            if l1 == l2
                ep = zeros(nz,1); ep(l1) = h;
                Fp = eval_rbc_F(ep,  i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                Fm = eval_rbc_F(-ep, i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                hess = (Fp + Fm - 2*F0(:,i)) / h^2;
            else
                epp = zeros(nz,1); epp(l1)=h; epp(l2)=h;
                epm = zeros(nz,1); epm(l1)=h; epm(l2)=-h;
                emp = zeros(nz,1); emp(l1)=-h; emp(l2)=h;
                emm = zeros(nz,1); emm(l1)=-h; emm(l2)=-h;
                Fpp = eval_rbc_F(epp, i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                Fpm = eval_rbc_F(epm, i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                Fmp = eval_rbc_F(emp, i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                Fmm = eval_rbc_F(emm, i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_n, gh_w);
                hess = (Fpp - Fpm - Fmp + Fmm) / (4*h^2);
            end

            rows = (i-1)*n*n_unique + (m_idx-1)*n + 1 : (i-1)*n*n_unique + m_idx*n;
            hvec(rows) = hess;
        end
    end
end
end


function [ri, vi, l1, l2] = linear_idx_to_Q(m, ns, n, nz)
% Map linear index m to (regime, variable, l1, l2) in Q.
n_unique = nz*(nz+1)/2;
per_regime = n * n_unique;

ri = ceil(m / per_regime);
rem1 = m - (ri-1)*per_regime;
vi = ceil(rem1 / n_unique);
rem2 = rem1 - (vi-1)*n_unique;

cnt = 0;
for ll1 = 1:nz
    for ll2 = ll1:nz
        cnt = cnt + 1;
        if cnt == rem2
            l1 = ll1; l2 = ll2;
            return;
        end
    end
end
end
