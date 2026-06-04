function [css, kss, zss] = rbc_steady_state(par, mu_bar)
% Compute steady state of the stationary RBC model.
% From Foerster et al. (2016) eq (14), page 645.

alpha = par.alpha; beta = par.beta; nu = par.nu; delta = par.delta;

zss = exp(mu_bar);

temp = exp((1-nu)*mu_bar)/beta - 1 + delta;
kss = ((1/alpha) * exp((alpha-1)*mu_bar) * temp)^(1/(alpha-1));

css = kss * ((1/alpha)*temp + 1 - delta - exp(mu_bar));
end
