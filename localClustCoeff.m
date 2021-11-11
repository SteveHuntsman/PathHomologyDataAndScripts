function y = localClustCoeff(A)

% Input: digraph adjacency matrix A
% Based on https://doi.org/10.1103/PhysRevE.76.026107 (uses eqn 8 of the
% version at https://arxiv.org/abs/physics/0612169)

if any(diag(A)), warning('killing diagonal'); A = A-diag(diag(A)); end

d_tot = sum(A,1)'+sum(A,2); % column
d_bi = diag(A^2);           % column
y = .5*diag((A+A')^3)./(d_tot.*(d_tot-1)-2*d_bi);