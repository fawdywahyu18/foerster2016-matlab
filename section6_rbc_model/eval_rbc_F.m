function Fval = eval_rbc_F(z, regime_i, sol, Q, par, P, mu_vec, rho_vec, sigma_vec, method, gh_nodes, gh_weights)
% Evaluate F_i(z) = sum_j p_{ij} * integral f(...) d mu(eps_{t+1})
% using Gauss-Hermite quadrature for the eps_{t+1} integration.
%
% z = [xhat_{t-1,1}; xhat_{t-1,2}; eps_t; chi]  (4x1)
% Q{regime, var}: nz x nz symmetric second-order coefficient matrices

ns = size(P,1); ny = 1; nx = 2; nz = 4;
yss = sol.yss; xss = sol.xss;
mu_bar = sol.mu_bar; rho_bar = sol.rho_bar; sigma_bar = sol.sigma_bar;

chi = z(4);

Ay_i = [sol.G{regime_i}, sol.Psi{regime_i}, sol.d{regime_i}];
Ax_i = [sol.H{regime_i}, sol.Phi{regime_i}, sol.cc{regime_i}];

y_t  = yss    + Ay_i * z    + 0.5 * z' * Q{regime_i,1} * z;
xt1  = xss(1) + Ax_i(1,:)*z + 0.5 * z' * Q{regime_i,2} * z;
xt2  = xss(2) + Ax_i(2,:)*z + 0.5 * z' * Q{regime_i,3} * z;
x_t  = [xt1; xt2];

xlag  = xss + z(1:2);
eps_t = z(3);

Fval = zeros(3,1);
for j = 1:ns
    if strcmp(method, 'partition')
        theta_j = [mu_bar + chi*(mu_vec(j)-mu_bar); rho_vec(j); sigma_vec(j)];
        theta_i = [mu_bar + chi*(mu_vec(regime_i)-mu_bar); rho_vec(regime_i); sigma_vec(regime_i)];
    else
        theta_j = [mu_bar + chi*(mu_vec(j)-mu_bar);
                    rho_bar + chi*(rho_vec(j)-rho_bar);
                    sigma_bar + chi*(sigma_vec(j)-sigma_bar)];
        theta_i = [mu_bar + chi*(mu_vec(regime_i)-mu_bar);
                    rho_bar + chi*(rho_vec(regime_i)-rho_bar);
                    sigma_bar + chi*(sigma_vec(regime_i)-sigma_bar)];
    end

    for k = 1:length(gh_nodes)
        eps_tp1 = gh_nodes(k);

        xhat_t  = x_t - xss;
        z_next  = [xhat_t; chi*eps_tp1; chi];

        Ay_j = [sol.G{j}, sol.Psi{j}, sol.d{j}];
        y_tp1 = yss + Ay_j * z_next + 0.5 * z_next' * Q{j,1} * z_next;

        fval = rbc_eval_f(y_tp1, y_t, x_t, xlag, chi*eps_tp1, eps_t, theta_j, theta_i, par);
        Fval = Fval + P(regime_i,j) * gh_weights(k) * fval;
    end
end
end
