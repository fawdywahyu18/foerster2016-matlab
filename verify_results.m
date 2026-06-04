%% Verification script: compare all computed results against Foerster et al. (2016)
clear; clc;
fprintf('=============================================================\n');
fprintf('  VERIFICATION: Foerster et al. (2016) Replication\n');
fprintf('=============================================================\n\n');

tol = 0.02;
n_pass = 0; n_fail = 0; n_total = 0;

%% ===================================================================
%  SECTION 5: Simple Inflation Model (Table 1)
%  ===================================================================
fprintf('--- SECTION 5: Simple Model ---\n\n');

P_simple = [0.95, 0.05; 0.15, 0.85];
c0 = [1.25; 0.96];
c1 = [0.0; 0.0];
gsigma = [0.1; 0.6];

pbar_s = ergodic_probs(P_simple);
c0bar = pbar_s' * c0;
sigmabar = pbar_s' * gsigma;

% Partition method: exact solution
xi_part = [-gsigma(1)/c0(1); -gsigma(2)/c0(2)];
check_val('Partition xi(1)', xi_part(1), -0.08, tol);
check_val('Partition xi(2)', xi_part(2), -0.625, tol);

% Naive method 1st order
xi_naive = -sigmabar / c0bar;
check_val('Naive 1st xi', xi_naive, -0.191083, tol);

% Euler equation errors
gpitm1 = 0.01; gepst = 1.0;

ee_part = 0;
for st = 1:2
    gpit = xi_part(st)*gepst;
    Ea = P_simple(st,:) * [0;0];
    ee = c0(st)*gpit + c1(st)*gpitm1 - Ea*gpit + gsigma(st)*gepst;
    ee_part = ee_part + abs(ee);
end
ee_part = ee_part / 2;
if ee_part < 1e-15
    ee_part_log = -Inf;
    fprintf('  PASS: Partition EE = -Inf (exact solution)\n');
    n_pass = n_pass + 1;
else
    ee_part_log = log10(ee_part);
    fprintf('  FAIL: Partition EE = %.4f (expected -Inf)\n', ee_part_log);
    n_fail = n_fail + 1;
end
n_total = n_total + 1;

ee_naive1 = 0;
xi_n1 = [xi_naive; xi_naive];
for st = 1:2
    gpit = xi_n1(st)*gepst;
    Ea = P_simple(st,:) * [0;0];
    ee = c0(st)*gpit + c1(st)*gpitm1 - Ea*gpit + gsigma(st)*gepst;
    ee_naive1 = ee_naive1 + abs(ee);
end
ee_naive1_log = log10(ee_naive1/2);
check_val('Naive 1st EE (log10)', ee_naive1_log, -0.5564, tol);

% Naive 2nd order
A_sys = c0bar * eye(2) - P_simple;
gxi_naive = -sigmabar / c0bar;
rhs_b5 = -((c0 - c0bar)*gxi_naive + (gsigma - sigmabar));
b5_half = A_sys \ rhs_b5;
b5 = 2 * b5_half;

ee_naive2 = 0;
gchi = 1.0;
for st = 1:2
    Ea = P_simple(st,:) * [0;0];
    gpit = xi_n1(st)*gepst + 0.5*b5(st)*gepst*gchi;
    ee = c0(st)*gpit + c1(st)*gpitm1 - Ea*gpit + gsigma(st)*gepst;
    ee_naive2 = ee_naive2 + abs(ee);
end
ee_naive2_log = log10(ee_naive2/2);
check_val('Naive 2nd EE (log10)', ee_naive2_log, -1.3691, tol);

%% ===================================================================
%  SECTION 6: RBC Model
%  ===================================================================
fprintf('\n--- SECTION 6: RBC Model ---\n\n');

par.alpha = 0.33; par.beta = 0.9976; par.nu = -1; par.delta = 0.025;
mu_vec = [0.0274; -0.0337]; rho_vec = [0.1; 0.0]; sigma_vec = [0.0072; 0.0216];
P = [0.75, 0.25; 0.50, 0.50];

pbar = ergodic_probs(P);
mu_bar = pbar' * mu_vec;

[css, kss, zss] = rbc_steady_state(par, mu_bar);
check_val('c_tilde_ss', css, 2.08259, 0.005);
check_val('k_tilde_ss', kss, 22.1504, 0.01);
check_val('z_tilde_ss', zss, 1.007, 0.001);

% First-order partition
fprintf('\n  Solving partition 1st order...\n');
sol_part = solve_rbc_first_order(par, P, mu_vec, rho_vec, sigma_vec, 'partition');

check_val('Part G1(1,1)', sol_part.G{1}(1), 0.0405, tol);
check_val('Part G1(1,2)', sol_part.G{1}(2), 0.1264, tol);
check_val('Part H1(1,1)', sol_part.H{1}(1,1), 0.9692, tol);
check_val('Part H1(1,2)', sol_part.H{1}(1,2), -2.1406, tol);
check_val('Part H1(2,1)', sol_part.H{1}(2,1), 0.0, tol);
check_val('Part H1(2,2)', sol_part.H{1}(2,2), 0.1, tol);
check_val('Part Psi1', sol_part.Psi{1}, 0.0091, tol);
check_val('Part Phi1(1)', sol_part.Phi{1}(1), -0.1552, tol);
check_val('Part Phi1(2)', sol_part.Phi{1}(2), 0.0072, tol);

check_val('Part G2(1,2)', sol_part.G{2}(2), 0.0, tol);
check_val('Part H2(1,2)', sol_part.H{2}(1,2), 0.0, tol);
check_val('Part H2(2,2)', sol_part.H{2}(2,2), 0.0, tol);
check_val('Part Psi2', sol_part.Psi{2}, 0.0268, tol);

% First-order naive
fprintf('\n  Solving naive 1st order...\n');
sol_naive = solve_rbc_first_order(par, P, mu_vec, rho_vec, sigma_vec, 'naive');

check_val('Naive G1(1,1)', sol_naive.G{1}(1), 0.0406, tol);
check_val('Naive G1(1,2)', sol_naive.G{1}(2), 0.0836, tol);
check_val('Naive H1(1,1)', sol_naive.H{1}(1,1), 0.9692, tol);
check_val('Naive H1(1,2)', sol_naive.H{1}(1,2), -1.4264, tol);
check_val('Naive H1(2,2)', sol_naive.H{1}(2,2), 0.0667, tol);
check_val('Naive Psi1', sol_naive.Psi{1}, 0.0152, tol);

% Verify naive coefficients are same across regimes
diff_G = norm(sol_naive.G{1} - sol_naive.G{2});
diff_H = norm(sol_naive.H{1} - sol_naive.H{2});
diff_Psi = norm(sol_naive.Psi{1} - sol_naive.Psi{2});
diff_Phi = norm(sol_naive.Phi{1} - sol_naive.Phi{2});
fprintf('  Naive regime-independence check (should be ~0):\n');
fprintf('    |G1-G2| = %.2e, |H1-H2| = %.2e, |Psi1-Psi2| = %.2e, |Phi1-Phi2| = %.2e\n', ...
    diff_G, diff_H, diff_Psi, diff_Phi);
if max([diff_G, diff_H, diff_Psi, diff_Phi]) < 1e-8
    fprintf('    PASS: Naive 1st-order is regime-independent\n');
    n_pass = n_pass + 1;
else
    fprintf('    FAIL: Naive 1st-order NOT regime-independent\n');
    n_fail = n_fail + 1;
end
n_total = n_total + 1;

%% Summary
fprintf('\n=============================================================\n');
fprintf('  VERIFICATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  Passed: %d / %d\n', n_pass, n_total);
fprintf('  Failed: %d / %d\n', n_fail, n_total);
if n_fail == 0
    fprintf('  ALL CHECKS PASSED!\n');
else
    fprintf('  Some checks failed - review output above.\n');
end
fprintf('=============================================================\n');


function check_val(name, computed, expected, tol)
    persistent np nf nt;
    if isempty(np), np = 0; nf = 0; nt = 0; end

    err = abs(computed - expected);
    rel_err = err / max(abs(expected), 1e-10);
    if rel_err < tol || err < tol
        fprintf('  PASS: %-25s = %10.4f  (expected %10.4f, err=%.1e)\n', name, computed, expected, err);
        np = np + 1;
    else
        fprintf('  FAIL: %-25s = %10.4f  (expected %10.4f, err=%.1e)\n', name, computed, expected, err);
        nf = nf + 1;
    end
    nt = nt + 1;

    assignin('caller', 'n_pass', evalin('caller','n_pass') + (rel_err < tol || err < tol));
    assignin('caller', 'n_fail', evalin('caller','n_fail') + ~(rel_err < tol || err < tol));
    assignin('caller', 'n_total', evalin('caller','n_total') + 1);
end
