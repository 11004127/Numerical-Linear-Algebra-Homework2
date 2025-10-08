

% Input matrix A
A = input('Entry the value of square matrix A: ');
disp('Matrix A =');
disp(A);

[m, n] = size(A);
if m < n
    error('Matrix must satisfy m ≥ n for QR decomposition!');
end

% Generate b
mode = input('Generate b automatically? (1=yes, 0=manual): ');
if mode == 1
    x_hat = randn(n,1); % ground-truth
    b = A * x_hat;       % corresponding b
    disp('Generated ground-truth x̂ ='); disp(x_hat);
else
    b = input('Enter the value of vector b: ');
end

% QR decomposition using Householder
Q = eye(m);
R = A;
for k = 1:n
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

    % Update R and Q
    R(k:end, k:end) = R(k:end, k:end) - 2 * v * (v' * R(k:end, k:end));
    Q(:, k:end) = Q(:, k:end) - 2 * Q(:, k:end) * (v * v');
end
R(abs(R) < 1e-12) = 0;

% Least squares using QR
y = Q' * b;
R1 = R(1:n, 1:n);
y1 = y(1:n);
x_QR = R1 \ y1;

disp('Least squares solution x_QR (using QR) =');
disp(x_QR);

% SVD decomposition
[U, S, V] = svd(A, 'econ');

% Manually compute S_inv
[nRows, nCols] = size(S);
S_inv = zeros(nCols, nRows);
for i = 1:min(nRows, nCols)
    if S(i,i) > 1e-12
        S_inv(i,i) = 1 / S(i,i);
    end
end

% Least squares using SVD
x_svd = V * S_inv * (U' * b);

disp('Least squares solution x_svd (using SVD) =');
disp(x_svd);

% Compare results
diff_norm = norm(x_svd - x_QR, 2);
fprintf('Difference between SVD and QR solution ||x_svd - x_QR||_2 = %.3e\n', diff_norm);

% Residual check
residual_QR = norm(A*x_QR - b, 2);
residual_svd = norm(A*x_svd - b, 2);
fprintf('Residual norm QR ||A*x_QR - b||_2 = %.3e\n', residual_QR);
fprintf('Residual norm SVD ||A*x_svd - b||_2 = %.3e\n', residual_svd);

% Optional: error to ground-truth
if exist('x_hat','var')
    abs_err_QR = norm(x_QR - x_hat, 2);
    abs_err_SVD = norm(x_svd - x_hat, 2);
    fprintf('Absolute error QR ||x_QR - x̂||_2 = %.3e\n', abs_err_QR);
    fprintf('Absolute error SVD ||x_svd - x̂||_2 = %.3e\n', abs_err_SVD);
end
