L = [3, 0, 0; 1, 2, 0; 4, 5, 6];
b = [6; 8; 24];
x = trilsolve(L, b);

expected_solution = L \ b;
tolerance = 1e-6;
if norm(x - expected_solution) < tolerance
    disp('Test passed: Solutions match within tolerance.');
else
    disp('Test failed: Solutions do not match within tolerance.');
end

function x = trilsolve(L, b)
    N = size(L, 1);
    x = zeros(N, 1);

    for n = 1:N
        % Define matrix and vector blocks.
        L11 = L(n, n);
        L21 = L((n+1):end, n);
        b1 = b((n+1):end);
        b2 = b(n);

        % solve x(i).
        x(n) = b2 / L11;

        % update vector b
        b((n+1):end) = b1 - L21 * x(n);
    end
end

