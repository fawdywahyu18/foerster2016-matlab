function [ee_mean, ee_log10] = compute_ee_rbc(sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, order)
% Compute Euler equation errors for the RBC model via simulation.
% Follows the procedure in Foerster et al. (2016), Section 6.3.
%
% Simulates 10,000 periods (burn-in 1,000), then for each of the
% remaining 9,000 draws, computes the EE residual using Monte Carlo
% integration over eps_{t+1} and s_{t+1}.

ns = size(P,1); ny = 1; nx = 2; nz = 4;

yss = sol.yss; xss = sol.xss;
css = sol.css; kss = sol.kss; zss = sol.zss;

alpha = par.alpha; beta = par.beta; nu = par.nu; delta = par.delta;

T_total = 10000;
T_burn  = 1000;
N_mc    = 10000;

if nargin < 9
    order = 1;
end

if isempty(Q)
    Q = cell(ns, 3);
    for i = 1:ns
        for v = 1:3
            Q{i,v} = zeros(nz);
        end
    end
end

A_y = cell(ns,1); A_x = cell(ns,1);
for i = 1:ns
    A_y{i} = [sol.G{i}, sol.Psi{i}, sol.d{i}];
    A_x{i} = [sol.H{i}, sol.Phi{i}, sol.cc{i}];
end

rng(42);
eps_sim = randn(T_total, 1);

cum_P = cumsum(P, 2);
s_sim = zeros(T_total, 1);
s_sim(1) = 1 + (rand < sol.pbar(2));
for t = 2:T_total
    u = rand;
    s_sim(t) = find(u <= cum_P(s_sim(t-1),:), 1);
end

khat = zeros(T_total, 1);
zhat = zeros(T_total, 1);
chat = zeros(T_total, 1);

for t = 1:T_total
    si = s_sim(t);
    if t == 1
        z_vec = [0; 0; eps_sim(t); 1];
    else
        z_vec = [khat(t-1); zhat(t-1); eps_sim(t); 1];
    end

    chat(t) = A_y{si} * z_vec;
    khat_t  = A_x{si}(1,:) * z_vec;
    zhat_t  = A_x{si}(2,:) * z_vec;

    if order >= 2
        chat(t) = chat(t) + 0.5 * z_vec' * Q{si,1} * z_vec;
        khat_t  = khat_t  + 0.5 * z_vec' * Q{si,2} * z_vec;
        zhat_t  = zhat_t  + 0.5 * z_vec' * Q{si,3} * z_vec;
    end

    khat(t) = khat_t;
    zhat(t) = zhat_t;
end

c_sim = css + chat;
k_sim = kss + khat;
z_sim = zss + zhat;

c_sim = max(c_sim, 1e-10);
k_sim = max(k_sim, 1e-10);
z_sim = max(z_sim, 1e-10);

eps_mc = randn(N_mc, 1);

ee_vals = zeros(T_total - T_burn, 1);

for tt = 1:(T_total - T_burn)
    t = T_burn + tt;
    si = s_sim(t);
    ct = c_sim(t);
    kt = k_sim(t);
    zt = z_sim(t);

    khat_t = khat(t);
    zhat_t = zhat(t);

    rhs_sum = 0;

    for j = 1:ns
        mu_j    = mu_vec(j);
        rho_j   = rho_vec(j);
        sigma_j = sigma_vec(j);

        z_next_base = [khat_t; zhat_t];

        c_tp1_mc = zeros(N_mc, 1);
        z_tp1_mc = zeros(N_mc, 1);

        for m = 1:N_mc
            ep = eps_mc(m);

            z_tilde_tp1 = exp(mu_j + rho_j*(log(zt) - mu_j) + sigma_j*ep);

            z_vec_next = [khat_t; zhat_t; ep; 1];
            c_hat_tp1 = A_y{j} * z_vec_next;
            if order >= 2
                c_hat_tp1 = c_hat_tp1 + 0.5 * z_vec_next' * Q{j,1} * z_vec_next;
            end

            c_tp1_mc(m) = max(css + c_hat_tp1, 1e-10);
            z_tp1_mc(m) = z_tilde_tp1;
        end

        integrand = c_tp1_mc.^(nu-1) .* ...
                    (alpha * z_tp1_mc.^(1-alpha) .* kt^(alpha-1) + 1 - delta);
        rhs_sum = rhs_sum + P(si,j) * mean(integrand);
    end

    ee_val = 1 - beta * zt^(nu-1) * rhs_sum / ct^(nu-1);
    ee_vals(tt) = abs(ee_val);
end

ee_mean = mean(ee_vals);
if ee_mean < 1e-20
    ee_log10 = -Inf;
else
    ee_log10 = log10(ee_mean);
end
end
