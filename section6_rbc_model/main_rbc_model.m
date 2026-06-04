%% Foerster et al. (2016) - Section 6: Markov-Switching RBC Model
% Replicates Table 3: Euler-equation errors (base-10 log absolute value)
%
% Model: 3 equilibrium equations (Euler, resource constraint, technology)
%   y_t = c_tilde_t,  x_t = [k_tilde_t, z_tilde_t]
%   n_y = 1, n_x = 2, n_e = 1, n_s = 2
%
% Two methods: partition perturbation vs naive perturbation
% Orders: 1st and 2nd

clear; clc;
fprintf('=============================================================\n');
fprintf('  Foerster et al. (2016) - Section 6: RBC Model\n');
fprintf('  Replication of Table 3\n');
fprintf('=============================================================\n\n');

%% Parameters (Table 2)
par.alpha = 0.33;
par.beta  = 0.9976;
par.nu    = -1;
par.delta = 0.025;

mu_vec    = [0.0274; -0.0337];
rho_vec   = [0.1; 0.0];
sigma_vec = [0.0072; 0.0216];

P = [0.75, 0.25;
     0.50, 0.50];

ns = size(P,1);

%% Ergodic probabilities and means
pbar = ergodic_probs(P);
mu_bar    = pbar' * mu_vec;
rho_bar   = pbar' * rho_vec;
sigma_bar = pbar' * sigma_vec;

fprintf('Ergodic probabilities: [%.4f, %.4f]\n', pbar(1), pbar(2));
fprintf('Ergodic means: mu_bar=%.6f, rho_bar=%.6f, sigma_bar=%.6f\n\n', ...
    mu_bar, rho_bar, sigma_bar);

%% Steady State
[css, kss, zss] = rbc_steady_state(par, mu_bar);
fprintf('Steady state:\n');
fprintf('  c_tilde_ss = %.5f  (paper: 2.08259)\n', css);
fprintf('  k_tilde_ss = %.4f  (paper: 22.1504)\n', kss);
fprintf('  z_tilde_ss = %.4f  (paper: 1.007)\n\n', zss);

%% ===================================================================
%  PARTITION PERTURBATION METHOD (1st order)
%  ===================================================================
fprintf('--- PARTITION PERTURBATION: 1st order ---\n');
tic;
sol_part = solve_rbc_first_order(par, P, mu_vec, rho_vec, sigma_vec, 'partition');
t_part1 = toc;

fprintf('  Slope residual: %.2e\n', sol_part.slope_residual);
fprintf('  MSS stable: %d (max eigenvalue: %.4f)\n', sol_part.is_stable, sol_part.max_eig);
fprintf('  Time: %.2f sec\n\n', t_part1);

fprintf('  Regime 1 first-order solution:\n');
A1 = [sol_part.G{1}, sol_part.Psi{1}, sol_part.d{1};
      sol_part.H{1}, sol_part.Phi{1}, sol_part.cc{1}];
fprintf('    [c_hat]   [%7.4f %7.4f %7.4f %7.4f] [k_hat_{t-1}]\n', A1(1,:));
fprintf('    [k_hat] = [%7.4f %7.4f %7.4f %7.4f] [z_hat_{t-1}]\n', A1(2,:));
fprintf('    [z_hat]   [%7.4f %7.4f %7.4f %7.4f] [eps_t; 1   ]\n\n', A1(3,:));

fprintf('  Paper regime 1:\n');
fprintf('    [0.0405  0.1264  0.0091  0.0000]\n');
fprintf('    [0.9692 -2.1406 -0.1552 -0.3720]\n');
fprintf('    [0.0000  0.1000  0.0072  0.0184]\n\n');

fprintf('  Regime 2 first-order solution:\n');
A2 = [sol_part.G{2}, sol_part.Psi{2}, sol_part.d{2};
      sol_part.H{2}, sol_part.Phi{2}, sol_part.cc{2}];
fprintf('    [c_hat]   [%7.4f %7.4f %7.4f %7.4f] [k_hat_{t-1}]\n', A2(1,:));
fprintf('    [k_hat] = [%7.4f %7.4f %7.4f %7.4f] [z_hat_{t-1}]\n', A2(2,:));
fprintf('    [z_hat]   [%7.4f %7.4f %7.4f %7.4f] [eps_t; 1   ]\n\n', A2(3,:));

fprintf('  Paper regime 2:\n');
fprintf('    [0.0405  0.0000  0.0268 -0.0968]\n');
fprintf('    [0.9692  0.0000 -0.4649  0.9227]\n');
fprintf('    [0.0000  0.0000  0.0217 -0.0410]\n\n');

%% ===================================================================
%  NAIVE PERTURBATION METHOD (1st order)
%  ===================================================================
fprintf('--- NAIVE PERTURBATION: 1st order ---\n');
tic;
sol_naive = solve_rbc_first_order(par, P, mu_vec, rho_vec, sigma_vec, 'naive');
t_naive1 = toc;

fprintf('  Slope residual: %.2e\n', sol_naive.slope_residual);
fprintf('  MSS stable: %d (max eigenvalue: %.4f)\n', sol_naive.is_stable, sol_naive.max_eig);
fprintf('  Time: %.2f sec\n\n', t_naive1);

fprintf('  First-order solution (same for both regimes):\n');
An = [sol_naive.G{1}, sol_naive.Psi{1}, sol_naive.d{1};
      sol_naive.H{1}, sol_naive.Phi{1}, sol_naive.cc{1}];
fprintf('    [c_hat]   [%7.4f %7.4f %7.4f %7.4f] [k_hat_{t-1}]\n', An(1,:));
fprintf('    [k_hat] = [%7.4f %7.4f %7.4f %7.4f] [z_hat_{t-1}]\n', An(2,:));
fprintf('    [z_hat]   [%7.4f %7.4f %7.4f %7.4f] [eps_t; 1   ]\n\n', An(3,:));

fprintf('  Paper values:\n');
fprintf('    [0.0406  0.0836  0.0152  0.0314]\n');
fprintf('    [0.9692 -1.4264 -0.2586 -0.4169]\n');
fprintf('    [0.0000  0.0667  0.0121  0.0191]\n\n');

%% ===================================================================
%  PARTITION PERTURBATION METHOD (2nd order)
%  ===================================================================
fprintf('--- PARTITION PERTURBATION: 2nd order ---\n');
fprintf('  Solving second-order system (this may take a moment)...\n');
tic;
Q_part = solve_rbc_second_order(sol_part, par, P, mu_vec, rho_vec, sigma_vec, 'partition');
t_part2 = toc;
fprintf('  Time: %.2f sec\n\n', t_part2);

display_second_order(Q_part, sol_part, 1, 'Partition regime 1');
display_second_order(Q_part, sol_part, 2, 'Partition regime 2');

%% ===================================================================
%  NAIVE PERTURBATION METHOD (2nd order)
%  ===================================================================
fprintf('--- NAIVE PERTURBATION: 2nd order ---\n');
fprintf('  Solving second-order system...\n');
tic;
Q_naive = solve_rbc_second_order(sol_naive, par, P, mu_vec, rho_vec, sigma_vec, 'naive');
t_naive2 = toc;
fprintf('  Time: %.2f sec\n\n', t_naive2);

display_second_order(Q_naive, sol_naive, 1, 'Naive regime 1');

%% ===================================================================
%  EULER EQUATION ERRORS
%  ===================================================================
fprintf('=============================================================\n');
fprintf('  Computing Euler equation errors (simulation-based)\n');
fprintf('  T=10000, burn-in=1000, MC draws=10000\n');
fprintf('=============================================================\n\n');

fprintf('  Partition 1st order...\n');
tic;
[~, ee_part1] = compute_ee_rbc(sol_part, [], par, P, mu_vec, rho_vec, sigma_vec, 'partition', 1);
fprintf('    EE (log10) = %.4f   (paper: -3.01)   Time: %.1f sec\n\n', ee_part1, toc);

fprintf('  Partition 2nd order...\n');
tic;
[~, ee_part2] = compute_ee_rbc(sol_part, Q_part, par, P, mu_vec, rho_vec, sigma_vec, 'partition', 2);
fprintf('    EE (log10) = %.4f   (paper: -3.59)   Time: %.1f sec\n\n', ee_part2, toc);

fprintf('  Naive 1st order...\n');
tic;
[~, ee_naive1] = compute_ee_rbc(sol_naive, [], par, P, mu_vec, rho_vec, sigma_vec, 'naive', 1);
fprintf('    EE (log10) = %.4f   (paper: -2.48)   Time: %.1f sec\n\n', ee_naive1, toc);

fprintf('  Naive 2nd order...\n');
tic;
[~, ee_naive2] = compute_ee_rbc(sol_naive, Q_naive, par, P, mu_vec, rho_vec, sigma_vec, 'naive', 2);
fprintf('    EE (log10) = %.4f   (paper: -3.07)   Time: %.1f sec\n\n', ee_naive2, toc);

%% ===================================================================
%  TABLE 3 COMPARISON
%  ===================================================================
fprintf('=============================================================\n');
fprintf('  TABLE 3: Euler-equation errors (base-10 log absolute value)\n');
fprintf('=============================================================\n');
fprintf('                        Computed    Paper\n');
fprintf('  Partition 1st order:   %6.2f     -3.01\n', ee_part1);
fprintf('  Partition 2nd order:   %6.2f     -3.59\n', ee_part2);
fprintf('  Naive 1st order:       %6.2f     -2.48\n', ee_naive1);
fprintf('  Naive 2nd order:       %6.2f     -3.07\n', ee_naive2);
fprintf('=============================================================\n');


%% ===================================================================
%  HELPER FUNCTION
%  ===================================================================
function display_second_order(Q, sol, regime, label)
    nz = 4;
    fprintf('  %s second-order Hessian (z kron z columns, 3 rows):\n', label);
    fprintf('  z = [khat, zhat, eps, chi], full 4x4 blocks:\n');

    for v = 1:3
        names = {'c_hat', 'k_hat', 'z_hat'};
        fprintf('    %s Hessian Q(%d,%d):\n', names{v}, regime, v);
        for r = 1:nz
            fprintf('      [');
            for c = 1:nz
                fprintf('%8.4f', Q{regime,v}(r,c));
            end
            fprintf(']\n');
        end
    end
    fprintf('\n');
end
