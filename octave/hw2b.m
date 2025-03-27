1;

% Setup: 
% A is diagonal with eigenvalues being a distribution along [a, 2] such that 2/a = desired condition number.
% Note: this definition of condition number holds since A is setup to be normal
% b = 1^T / sqrt(N), N dim of A (sqrt(N) is included to normalize b)
% analytic solution => x_i = 1 / (sqrt(N) * lambda_i)
% x0 = 0 so r0 = b
% Create a vector of 50 elements equally spaced between 1 and 2

% Given:
% N (dim)
% condition number
% (optional) eigenvalue distribution
% NOTE: I think uniform is fine

function errors_normalized = cg_bounds(N, k)
  a = 2 / k; % get lhs bound on eigen
  rho = (sqrt(k) - 1) / (sqrt(k) + 1); % upper bound we're testing
  sigma = linspace(a, 2, N); %eigenvalues
  A = diag(sigma); % construct A
  b = ones(N,1) ./sqrt(N); % b chosen so x_true is as below
  x_true = ones(N,1) ./ (sqrt(N) * diag(A));
  errors_normalized = x_true; % TEMP so function type is accurate
  [xm, ~, ~, m, resvec] = pcg(A, b, 1e-16);
  e0 = sqrt(x_true' * A * x_true); % compute || e0 ||_A
  em = sqrt((x_true-xm)' * A * (x_true-xm)); % compute || em ||_A
  disp("\n\n\n");
  disp("em");
  disp(em);
  disp("\n\n");
  disp("e0");
  disp(e0);
  disp("\n\n");
  disp(em / e0);
  disp(rho ^ m);
  disp((em / e0) < rho^m)
end

cg_bounds(30, 2);
% Do:
% 0. Define b, x0.
% 1. compute rho
% 2. construct A with cond. num and eigenvalue distribution
% 3. Compute x_true
% 4. Do CG, get error at each iteration
% 5. Normalize error by e0
% 6. Plot rho^m (upper bound)
% 7. Plot Normalized errors

% 8. Repeat for 5 different condition numbers

% Create the diagonal matrix
A = diag(sigma);

% Display the size and a few elements to verify
size(A)
A(1,1)  % First diagonal element (should be 1)
A(25,25) % Middle diagonal element (should be around 1.5)
A(50,50) % Last diagonal element (should be 2)

cond(A)

% TODO: use paper to design A
