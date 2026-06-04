function D = rbc_jacobians(yss, xss, theta_j, theta_i, par)
% Compute numerical Jacobians of f at the steady-state evaluation point.
% Central differences with step size h.

h = 1e-7;
ny = length(yss); nx = length(xss); n = ny + nx; nth = 3;

base = @(yp,y,x,xl,ep,e,tp,ti) rbc_eval_f(yp,y,x,xl,ep,e,tp,ti,par);

D.Dyp   = numjac(@(v) base(v,yss,xss,xss,0,0,theta_j,theta_i), yss, n, ny, h);
D.Dy    = numjac(@(v) base(yss,v,xss,xss,0,0,theta_j,theta_i), yss, n, ny, h);
D.Dx    = numjac(@(v) base(yss,yss,v,xss,0,0,theta_j,theta_i), xss, n, nx, h);
D.Dxlag = numjac(@(v) base(yss,yss,xss,v,0,0,theta_j,theta_i), xss, n, nx, h);
D.Depsp = numjac(@(v) base(yss,yss,xss,xss,v,0,theta_j,theta_i), 0, n, 1, h);
D.Deps  = numjac(@(v) base(yss,yss,xss,xss,0,v,theta_j,theta_i), 0, n, 1, h);
D.Dthetap = numjac(@(v) base(yss,yss,xss,xss,0,0,v,theta_i), theta_j, n, nth, h);
D.Dtheta  = numjac(@(v) base(yss,yss,xss,xss,0,0,theta_j,v), theta_i, n, nth, h);
end

function J = numjac(f, x0, nf, nv, h)
J = zeros(nf, nv);
for k = 1:nv
    pert = zeros(size(x0)); pert(k) = h;
    J(:,k) = (f(x0+pert) - f(x0-pert)) / (2*h);
end
end
