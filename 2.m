

% Entry the value of square matrix A
A = input('Entry the value of square matrix A: ');
disp('Matrix A =');
disp(A);

% Settings
[m, n] = size(A);
if m < n
  error('Matrix must satisfy m ≥ n for QR decomposition!');
end

% Generate b
mode = input('Generate b automatically? (1=yes, 0=manual): ');
if mode == 1
  x_hat = randn(n,1); % real solution
  b = A * x_hat; % corresponding b
  disp('Generated ground-truth x̂ ='); disp(x_hat);
else
  b = input('Enter the value of vector b: ');
end

% QR decomposition
Q = eye(m);
R = A;
for k = 1
  x = R(k:end, k);
  if all(x == 0)
    continue;
  end
  sigma = norm(x);
  alpha = -sign(x(1));
  if alpha == 0
    alpha = -1;
  end
  v = x;
  v(1) = v(1) - alpha * sigma;
  v = v / norm(v);
  R(k:end, k:end) = R(k:end, k:end) - 2 * v * (v' * R(k:end, k:end));
  Q(:, k:end) = Q(:, k:end) - 2 * Q(:, k:end) * (v * v');
end
R(abs(R) < 1e-12) = 0;

% Least squares problem
y = Q' * b;
R1 = R(1:n, 1:n);  % Upper triangular matrix
y1 = y(1:n);
x_QR = R1 \ y1;    % Solve upper triangular equations

disp('Least squares solution x_QR =');
disp(x_QR);

% Validation
if exist('x_hat','var')
  abs_err = norm(x_QR - x_hat, 2);
  fprintf('Absolute error ||x_QR - x̂||₂ = %.3e\n', abs_err);
end

residual = norm(A * x_QR - b, 2);
fprintf('Residual norm ||A*x_QR - b||₂ = %.3e\n', residual);

fprintf('Factorization check ||A - QR||_F = %.3e\n', norm(A - Q*R, 'fro'));
fprintf('Orthogonality check ||Q''*Q - I||_F = %.3e\n', norm(Q'*Q - eye(m), 'fro'));
