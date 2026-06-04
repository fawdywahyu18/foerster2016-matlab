function fval = rbc_eval_f(yp, y, x, xlag, epsp, eps_val, theta_p, theta, par)
% Evaluate the 3 equilibrium conditions of the Foerster et al. (2016) RBC model.
%
% f1: Euler equation
% f2: Resource constraint
% f3: Technology process
%
% yp = [cp], y = [c], x = [k,z], xlag = [klag,zlag]
% epsp = eps_{t+1}, eps_val = eps_t
% theta_p = [mu_j, rho_j, sigma_j]  (future regime s_{t+1})
% theta   = [mu_i, rho_i, sigma_i]  (current regime s_t)

cp = yp(1); c = y(1);
k = x(1);   z = x(2);
klag = xlag(1); zlag = xlag(2);

alpha = par.alpha; beta = par.beta; nu = par.nu; delta = par.delta;

mu_j = theta_p(1); rho_j = theta_p(2); sigma_j = theta_p(3);
mu_i = theta(1);   rho_i = theta(2);   sigma_i = theta(3);

zp = exp(mu_j + rho_j*(log(z) - mu_j) + sigma_j*epsp);

fval = zeros(3,1);
fval(1) = -c^(nu-1) + beta * z^(nu-1) * cp^(nu-1) * ...
          (alpha * zp^(1-alpha) * k^(alpha-1) + 1 - delta);
fval(2) = -c - k*z + z^(1-alpha) * klag^alpha + (1-delta)*klag;
fval(3) = -log(z) + mu_i + rho_i*(log(zlag) - mu_i) + sigma_i*eps_val;
end
