n = 4;
F = randn(n, n);

A = F * F';

L = rchol(A);

disp(L)

tolerance = 1e-6;

is_equal = norm(A - L * L', 'fro') <= tolerance;

if is_equal
    disp('Cholesky decomposition successful.');
else
    disp('Cholesky decomposition failed.');
end

function L = rchol(A)
    [n, n] = size(A);
    
    if n == 1
        L = sqrt(A);
    else
        % Split A
        a11 = A(1, 1);
        a21 = A(2:end, 1);
        a12 = A(1, 2:end);
        A22 = A(2:end, 2:end);
        
        % Compute L2
        L2 = rchol(A22 - (a21 * a12) / a11);

        L = eye(n);
        L(2:end, 2:end) = L2;
        L(2:end, 1) = a21 / sqrt(a11);
        L(1, 1) = sqrt(a11);
    end
end
