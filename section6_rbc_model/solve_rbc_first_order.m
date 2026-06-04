function sol = solve_rbc_first_order(par, P, mu_vec, rho_vec, sigma_vec, method)
% Solve the first-order MSDSGE system for the RBC model.
%
% method: 'partition' or 'naive'
%
% Returns sol struct with:
%   G{i}   (ny x nx)  - slope of control on states
%   H{i}   (nx x nx)  - slope of states on lagged states
%   Psi{i} (ny x ne)  - impact of shock on control
%   Phi{i} (nx x ne)  - impact of shock on states
%   d{i}   (ny x 1)   - chi coefficient for control
%   cc{i}  (nx x 1)   - chi coefficient for states

ns = size(P,1); ny = 1; nx = 2; ne = 1; nth = 3;
n = ny + nx;

pbar = ergodic_probs(P);
mu_bar = pbar' * mu_vec;
rho_bar = pbar' * rho_vec;
sigma_bar = pbar' * sigma_vec;

[css, kss, zss] = rbc_steady_state(par, mu_bar);
yss = css;
xss = [kss; zss];

% --- Compute Jacobians for each (i,j) regime pair ---
Jac = cell(ns, ns);
for i = 1:ns
    for j = 1:ns
        if strcmp(method, 'partition')
            theta_j = [mu_bar; rho_vec(j); sigma_vec(j)];
            theta_i = [mu_bar; rho_vec(i); sigma_vec(i)];
        else
            theta_j = [mu_bar; rho_bar; sigma_bar];
            theta_i = [mu_bar; rho_bar; sigma_bar];
        end
        Jac{i,j} = rbc_jacobians(yss, xss, theta_j, theta_i, par);
    end
end

% ====================================================================
%  1. SLOPE COEFFICIENTS (quadratic system via fsolve)
% ====================================================================
%  For regime i:
%    sum_j p_{ij} * [Dxlag + Dyp*G_j*H_i + Dy*G_i + Dx*H_i] = 0
%
%  Unknowns: phi = [vec(S_1); vec(S_2)] where S_i = [G_i; H_i]

fun = @(phi) slope_residual(phi, Jac, P, ns, ny, nx);

phi0 = zeros(ns*n*nx, 1);
for i = 1:ns
    idx = (i-1)*n*nx;
    phi0(idx+ny*nx+1) = 1 - par.delta;
    if strcmp(method, 'partition')
        phi0(idx+ny*nx+4) = rho_vec(i);
    else
        phi0(idx+ny*nx+4) = rho_bar;
    end
end

opts = optimoptions('fsolve','Display','off','TolFun',1e-14,...
       'TolX',1e-14,'MaxFunEvals',100000,'MaxIter',10000,...
       'Algorithm','trust-region-dogleg');

[phi_sol, fval, exitflag] = fsolve(fun, phi0, opts);

if exitflag <= 0 || norm(fval) > 1e-8
    best_phi = phi_sol; best_norm = norm(fval);
    for trial = 1:100
        rng(trial);
        phi_try = phi0 + 0.05*randn(size(phi0));
        [phi_trial, fval_trial, ef] = fsolve(fun, phi_try, opts);
        if norm(fval_trial) < best_norm
            best_phi = phi_trial; best_norm = norm(fval_trial);
        end
        if best_norm < 1e-12, break; end
    end
    phi_sol = best_phi;
end

G = cell(ns,1); H = cell(ns,1);
for i = 1:ns
    idx = (i-1)*n*nx;
    S = reshape(phi_sol(idx+1:idx+n*nx), n, nx);
    G{i} = S(1:ny,:);
    H{i} = S(ny+1:end,:);
end

% Check MSS and keep only stable solution
[is_stable, max_eig] = check_mss(H, P);
if ~is_stable
    warning('First-order solution is NOT MSS stable (max eigenvalue = %.4f)', max_eig);
end

% ====================================================================
%  2. IMPACT COEFFICIENTS (linear system, decoupled by regime)
% ====================================================================
%  For regime i:
%    sum_j p_{ij} * [Dy*Psi_i + (Dyp*G_j + Dx)*Phi_i + Deps] = 0

Psi = cell(ns,1); Phi = cell(ns,1);
for i = 1:ns
    A_lhs = zeros(n, ny+nx);
    b_rhs = zeros(n, ne);
    for j = 1:ns
        A_lhs(:,1:ny)     = A_lhs(:,1:ny)     + P(i,j)*Jac{i,j}.Dy;
        A_lhs(:,ny+1:end) = A_lhs(:,ny+1:end) + P(i,j)*(Jac{i,j}.Dyp*G{j} + Jac{i,j}.Dx);
        b_rhs = b_rhs - P(i,j)*Jac{i,j}.Deps;
    end
    impact = A_lhs \ b_rhs;
    Psi{i} = impact(1:ny,:);
    Phi{i} = impact(ny+1:end,:);
end

% ====================================================================
%  3. CHI COEFFICIENTS (linear system, coupled across regimes)
% ====================================================================
%  For regime i:
%    sum_j p_{ij} * [Dy*d_i + (Dyp*G_j+Dx)*c_i + Dyp*d_j
%                    + Dthetap*dtheta_j + Dtheta*dtheta_i] = 0

AA = zeros(ns*n, ns*(ny+nx));
bb = zeros(ns*n, 1);

for i = 1:ns
    rows = (i-1)*n+1:i*n;

    if strcmp(method, 'partition')
        dtheta_i = [mu_vec(i) - mu_bar; 0; 0];
    else
        dtheta_i = [mu_vec(i)-mu_bar; rho_vec(i)-rho_bar; sigma_vec(i)-sigma_bar];
    end

    for j = 1:ns
        if strcmp(method, 'partition')
            dtheta_j = [mu_vec(j) - mu_bar; 0; 0];
        else
            dtheta_j = [mu_vec(j)-mu_bar; rho_vec(j)-rho_bar; sigma_vec(j)-sigma_bar];
        end

        cols_di = (i-1)*(ny+nx)+1:(i-1)*(ny+nx)+ny;
        cols_ci = (i-1)*(ny+nx)+ny+1:i*(ny+nx);
        cols_dj = (j-1)*(ny+nx)+1:(j-1)*(ny+nx)+ny;

        AA(rows, cols_di) = AA(rows, cols_di) + P(i,j)*Jac{i,j}.Dy;
        AA(rows, cols_ci) = AA(rows, cols_ci) + P(i,j)*(Jac{i,j}.Dyp*G{j} + Jac{i,j}.Dx);
        AA(rows, cols_dj) = AA(rows, cols_dj) + P(i,j)*Jac{i,j}.Dyp;

        bb(rows) = bb(rows) - P(i,j)*(Jac{i,j}.Dthetap*dtheta_j + Jac{i,j}.Dtheta*dtheta_i);
    end
end

chi_coeff = AA \ bb;
d = cell(ns,1); cc = cell(ns,1);
for i = 1:ns
    idx = (i-1)*(ny+nx);
    d{i}  = chi_coeff(idx+1:idx+ny);
    cc{i} = chi_coeff(idx+ny+1:idx+ny+nx);
end

% --- Store results ---
sol.G = G; sol.H = H;
sol.Psi = Psi; sol.Phi = Phi;
sol.d = d; sol.cc = cc;
sol.yss = yss; sol.xss = xss;
sol.css = css; sol.kss = kss; sol.zss = zss;
sol.mu_bar = mu_bar; sol.rho_bar = rho_bar; sol.sigma_bar = sigma_bar;
sol.pbar = pbar;
sol.Jac = Jac;
sol.is_stable = is_stable;
sol.max_eig = max_eig;
sol.slope_residual = norm(fval);
end


function res = slope_residual(phi, Jac, P, ns, ny, nx)
n = ny + nx;
G = cell(ns,1); H = cell(ns,1);
for i = 1:ns
    idx = (i-1)*n*nx;
    S = reshape(phi(idx+1:idx+n*nx), n, nx);
    G{i} = S(1:ny,:);
    H{i} = S(ny+1:end,:);
end

res = zeros(ns*n*nx, 1);
for i = 1:ns
    R = zeros(n, nx);
    for j = 1:ns
        R = R + P(i,j)*(Jac{i,j}.Dxlag + Jac{i,j}.Dyp*G{j}*H{i} ...
                        + Jac{i,j}.Dy*G{i} + Jac{i,j}.Dx*H{i});
    end
    res((i-1)*n*nx+1:i*n*nx) = R(:);
end
end
