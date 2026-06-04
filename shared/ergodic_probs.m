function pbar = ergodic_probs(P)
% Compute ergodic probabilities from transition matrix P.
% P(i,j) = Prob(s_{t+1}=j | s_t=i), rows sum to 1.
% Returns column vector of ergodic probabilities.

ns = size(P,1);
A = [P' - eye(ns); ones(1,ns)];
b = [zeros(ns,1); 1];
pbar = A \ b;
end
