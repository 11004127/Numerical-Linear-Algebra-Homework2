

% Entry the value of square matrix A
A = input('Entry the value of square matrix A: ');
disp('Matrix A =');
disp(A);

% Settings
[m, n] = size(A);
if m < n
  error('Matrix must satisfy m ≥ n for QR decomposition!');
end

Q = eye(m);
R = A;

for k = 1:n
    x = R(k:m, k); % Delete subvector
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
    v = v / norm(v); % Normalized reflection vector

    % Update the k:m columns and k:n columns of R (left-multiply by H)
    R(k:m, k:n) = R(k:m, k:n) - 2 * v * (v' * R(k:m, k:n));

    % Update the example block of Q (without building the entire H)
    Q(:, k:m) = Q(:, k:m) - 2 * (Q(:, k:m) * v) * v';
end

% Clear extremely small values
R(abs(R) < 1e-12) = 0;

disp('Q ='); disp(Q);
disp('R ='); disp(R);
disp('Q*R ='); disp(Q * R);

fprintf('||A - Q*R||_F = %.3e\n', norm(A - Q*R, 'fro'));
fprintf('||Q''*Q - I||_F = %.3e\n', norm(Q'*Q - eye(m), 'fro'));
fprintf('||tril(R,-1)||_F = %.3e\n', norm(tril(R,-1), 'fro'));
