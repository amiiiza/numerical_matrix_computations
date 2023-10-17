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
    n = length(L);

    if n == 1
        x = b / L;
    else
        % Define matrix and vector blocks.
        L11 = L(1, 1);
        L21 = L(2:end, 1);
        L22 = L(2:end, 2:end);

        b1 = b(1);
        b2 = b(2:end);
        % solve x1.
        x1 = b1 / L11;

        % solve x2 using recursive function call.
        S = L22;
        x2 = trilsolve(S, b2- L21 * x1);
        x = [x1; x2];
    end
end