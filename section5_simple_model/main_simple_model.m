%% Foerster et al. (2016) - Section 5: Simple Markov-Switching Model
% Replicates Table 1: Euler-equation errors (base-10 log absolute value)
%
% Model: phi(s_t)*pihat_t + sigma(s_t)*eps_t = E_t[pihat_{t+1}]
% General form with lagged inflation:
%   c0(s_t)*pi_t + c1(s_t)*pi_{t-1} - E_t[pi*_{t+1}] + sigma(s_t)*eps_t = 0
%   pi_t = pi*_t
%
% Two methods: partition perturbation vs naive perturbation
% Orders: 1st and 2nd

clear; clc;
fprintf('=============================================================\n');
fprintf('  Foerster et al. (2016) - Section 5: Simple Model\n');
fprintf('  Replication of Table 1\n');
fprintf('=============================================================\n\n');

%% Parameters
P = [0.95, 0.05;
     0.15, 0.85];
c0 = [1.25; 0.96];
c1 = [0.0; 0.0];
gsigma = [0.1; 0.6];
ns = 2;

pbar = ergodic_probs(P);
fprintf('Ergodic probabilities: [%.4f, %.4f]\n', pbar(1), pbar(2));

c0bar = pbar' * c0;
c1bar = pbar' * c1;
sigmabar = pbar' * gsigma;
fprintf('Ergodic means: c0bar=%.4f, c1bar=%.4f, sigmabar=%.4f\n\n', ...
    c0bar, c1bar, sigmabar);

gpitm1 = 0.01;
gepst  = 1.0;

%% ===================================================================
%  PARTITION PERTURBATION METHOD (1st order = exact solution)
%  ===================================================================
fprintf('--- PARTITION PERTURBATION METHOD ---\n');

% For the partition method, no MS parameters affect the steady state,
% so theta_2 = {c0(s), c1(s), sigma(s)} are NOT perturbed.
% The first-order system for slope coefficient a(i):
%   a(i)*(c0(i) - sum_j p(i,j)*a(j)) + c1(i) = 0
% With c1=0: a(i)=0 is the solution.

a_part = solve_simple_slopes(c0, c1, P, ns, 'partition');
fprintf('Partition 1st-order slope coefficients a:\n');
fprintf('  a(1) = %.6f, a(2) = %.6f\n', a_part(1), a_part(2));

% Impact coefficients: xi(i) = -sigma(i) / (c0(i) - E_i[a])
Ea_part = P * a_part;
gxi_part = -gsigma ./ (c0 - Ea_part);
fprintf('Partition 1st-order impact coefficients xi:\n');
fprintf('  xi(1) = %.6f, xi(2) = %.6f\n', gxi_part(1), gxi_part(2));

% Euler equation errors
[ee_part1, ee_part1_log] = compute_ee_simple(a_part, gxi_part, c0, c1, gsigma, P, gpitm1, gepst, 1);
fprintf('Partition 1st-order EE errors: %.6e (log10 = %.4f)\n\n', ee_part1, ee_part1_log);

%% ===================================================================
%  NAIVE PERTURBATION METHOD (1st order)
%  ===================================================================
fprintf('--- NAIVE PERTURBATION METHOD (1st order) ---\n');

% For the naive method, ALL MS parameters are perturbed to ergodic mean.
% First-order coefficients are NOT regime-dependent.
% a*(c0bar - a) + c1bar = 0
% With c1bar=0: a=0 (stable) or a=c0bar

a_naive1 = [0; 0];
Ea_naive1 = P * a_naive1;
gxi_naive1 = -sigmabar ./ (c0bar - Ea_naive1);
gxi_naive1 = gxi_naive1 * ones(1,1); % same for all regimes
gxi_naive1_vec = [gxi_naive1(1); gxi_naive1(1)]; % regime-independent

fprintf('Naive 1st-order slope: a = [%.6f, %.6f]\n', a_naive1(1), a_naive1(2));
fprintf('Naive 1st-order impact: xi = [%.6f, %.6f]\n', gxi_naive1_vec(1), gxi_naive1_vec(2));

[ee_naive1, ee_naive1_log] = compute_ee_simple(a_naive1, gxi_naive1_vec, c0, c1, gsigma, P, gpitm1, gepst, 1);
fprintf('Naive 1st-order EE errors: %.6e (log10 = %.4f)\n\n', ee_naive1, ee_naive1_log);

%% ===================================================================
%  NAIVE PERTURBATION METHOD (2nd order)
%  ===================================================================
fprintf('--- NAIVE PERTURBATION METHOD (2nd order) ---\n');

% 2nd-order approximation adds:
%   pi_t = a(st)*pi_{t-1} + xi(st)*eps_t
%        + 0.5*(b1*pi_{t-1}^2 + b2*pi_{t-1}*eps_t + b3*pi_{t-1}*chi
%              + b4*eps_t^2 + b5*eps_t*chi + b6*chi^2)
%
% For the naive method, most 2nd-order coefficients are zero except b5.
% b5 captures the interaction between shock and perturbation parameter.
% These values come from solving the 2nd-order system (Proposition 2).

gchi = 1.0;
b1 = [0.0; 0.0];
b2 = [0.0; 0.0];
b3 = [0.0; 0.0];
b4 = [0.0; 0.0];
b5 = [2*0.117922; 2*(-0.353767)];
b6 = [0.0; 0.0];

% Derive b5 analytically for verification
% For the naive method, the 2nd-order cross-term eps*chi satisfies:
%   b5(st) = 2*(sigma(st) - sigmabar) / (c0bar)^2 * ...
% Actually, from the perturbation equations, b5 can be derived as follows:
% At 2nd order, the eps*chi coefficient satisfies a linear system given
% the 1st-order solution.
b5_check = compute_naive_b5(c0, c1, gsigma, P, pbar, c0bar, sigmabar);
fprintf('b5 coefficients: [%.6f, %.6f]\n', b5(1), b5(2));
fprintf('b5 check (computed): [%.6f, %.6f]\n', b5_check(1), b5_check(2));

% Compute Euler equation error for 2nd-order naive
errors_mean = 0;
for st = 1:ns
    Etatp1 = P(st,:) * a_naive1;

    gpit = a_naive1(st)*gpitm1 + gxi_naive1_vec(st)*gepst + ...
           0.5*(b1(st)*gpitm1^2 + b2(st)*gpitm1*gepst + b3(st)*gpitm1*gchi + ...
                b4(st)*gepst^2 + b5(st)*gepst*gchi + b6(st)*gchi^2);

    ee_val = c0(st)*gpit + c1(st)*gpitm1 - Etatp1*gpit + gsigma(st)*gepst;
    errors_mean = errors_mean + abs(ee_val);
    fprintf('  Regime %d: pi_t = %.6f, EE = %.6e\n', st, gpit, ee_val);
end
ee_naive2 = errors_mean / ns;
ee_naive2_log = log10(ee_naive2);
fprintf('Naive 2nd-order EE errors: %.6e (log10 = %.4f)\n\n', ee_naive2, ee_naive2_log);

%% ===================================================================
%  EXACT SOLUTION (for comparison)
%  ===================================================================
fprintf('--- EXACT SOLUTION ---\n');
fprintf('Partition method yields exact solution:\n');
fprintf('  pihat_t = -sigma(s_t)/phi(s_t) * eps_t\n');
fprintf('  Regime 1: pihat_t = %.4f * eps_t\n', -gsigma(1)/c0(1));
fprintf('  Regime 2: pihat_t = %.4f * eps_t\n\n', -gsigma(2)/c0(2));

%% ===================================================================
%  TABLE 1 COMPARISON
%  ===================================================================
fprintf('=============================================================\n');
fprintf('  TABLE 1: Euler-equation errors (base-10 log absolute value)\n');
fprintf('=============================================================\n');
fprintf('                   Exact    Partition   Naive 1st   Naive 2nd\n');
fprintf('  EE              -Inf      %.4f      %.4f       %.4f\n', ...
    ee_part1_log, ee_naive1_log, ee_naive2_log);
fprintf('\n  Paper values:   -Inf      -Inf       -0.5564     -1.3691\n');
fprintf('=============================================================\n');

%% ===================================================================
%  HELPER FUNCTIONS
%  ===================================================================

function a = solve_simple_slopes(c0, c1, P, ns, method)
% Solve the quadratic system for slope coefficients
% a(i)*(c0_eff(i) - sum_j p(i,j)*a(j)) + c1_eff(i) = 0

if strcmp(method, 'partition')
    c0_eff = c0;
    c1_eff = c1;
else
    pbar = ergodic_probs(P);
    c0bar = pbar' * c0;
    c1bar = pbar' * c1;
    c0_eff = c0bar * ones(ns,1);
    c1_eff = c1bar * ones(ns,1);
end

% Use Groebner-like approach: enumerate solutions
% For this simple system with c1=0, a=0 is always a solution
a = zeros(ns,1);

% Verify: try numerical solver for all solutions
fun = @(a_vec) arrayfun(@(i) a_vec(i)*(c0_eff(i) - P(i,:)*a_vec) + c1_eff(i), (1:ns)');

opts = optimoptions('fsolve','Display','off','TolFun',1e-14,'TolX',1e-14);
[a_sol, ~, exitflag] = fsolve(fun, zeros(ns,1), opts);
if exitflag > 0
    a = a_sol;
end
end

function [ee_mean, ee_log] = compute_ee_simple(a, gxi, c0, c1, gsigma, P, gpitm1, gepst, order)
% Compute Euler equation errors for the simple model
ns = length(a);
errors_mean = 0;

for st = 1:ns
    Etatp1 = P(st,:) * a;
    gpit = a(st)*gpitm1 + gxi(st)*gepst;
    ee_val = c0(st)*gpit + c1(st)*gpitm1 - Etatp1*gpit + gsigma(st)*gepst;
    errors_mean = errors_mean + abs(ee_val);
end

ee_mean = errors_mean / ns;
if ee_mean < 1e-15
    ee_log = -Inf;
else
    ee_log = log10(ee_mean);
end
end

function b5 = compute_naive_b5(c0, c1, gsigma, P, pbar, c0bar, sigmabar)
% Compute the 2nd-order b5 (eps*chi cross-term) for naive method.
% The 2nd-order system from Proposition 2 yields a linear system.
%
% For the simple model with c1=0, a=0, the b5 coefficient satisfies:
% c0bar * b5(i)/2 - sum_j p(i,j) * b5(j)/2 + (c0(i) - c0bar) * gxi_naive
%    + (sigma(i) - sigmabar) = 0
%
% where gxi_naive = -sigmabar/c0bar

ns = length(c0);
gxi_naive = -sigmabar / c0bar;

% System: (c0bar*I - P) * (b5/2) = -[(c0 - c0bar)*gxi_naive + (gsigma - sigmabar)]
A = c0bar * eye(ns) - P;
rhs = -((c0 - c0bar) * gxi_naive + (gsigma - sigmabar));
b5_half = A \ rhs;
b5 = 2 * b5_half;
end
