1;

% TODO: Create function which does cg_bounds, but adjusts definition of x_true for M = diag(A). It should also call pcg with this preconditioner

% TODO: Do as above but for tridiagonal (or triangular)

function [en, rho_m] = cg_bounds(N, k)
  a = 2 / k; % get lhs bound on eigen
  rho = (sqrt(k) - 1) / (sqrt(k) + 1); % upper bound we're testing
  sigma = linspace(a, 2, N); %eigenvalues
  A = diag(sigma); % construct A
  b = ones(N,1) ./sqrt(N); % b chosen so x_true is as below
  x_true = ones(N,1) ./ (sqrt(N) * diag(A));
  [xm, ~, ~, m, resvec] = pcg(A, b, 1e-14); % compute xm and store m (num iterations)
  e0 = sqrt(x_true' * A * x_true); % compute || e0 ||_A
  em = sqrt((x_true-xm)' * A * (x_true-xm)); % compute || em ||_A
  en = em / e0;
  rho_m = rho^m;
end

for i = [2 10 50 100 1000]
  [em, rm] = cg_bounds(30, i);
  disp(['condition number: ', num2str(i), ' em / e0 = ', num2str(em), ', rho^m = ', num2str(rm)]);
  disp(['em / e0 <= rho^m: ', mat2str(em <= rm)]);
end
