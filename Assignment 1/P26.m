% Test the LU decomposition function
A = [4, 3, 6; 6, 3, 12; 12, 6, 27];
[P, L, U] = lu(A);

% Verify the correctness of the decomposition
disp("Original A:");
disp(A);
disp("P:");
disp(P);
disp("L:");
disp(L);
disp("U:");
disp(U);
disp("P * L * U:");
disp(P * L * U);

function [P, L, U] = lu(A)
    [n, ~] = size(A);
    P = eye(n);
    L = eye(n);
    U = A;
    
    % Base case: A is 1x1 or empty
    if n <= 1
        return;
    end
    
    % Find the pivot element and exchange rows if necessary
    [~, pivot_row] = max(abs(U(:, 1)));
    if pivot_row ~= 1
        % Swap rows in P, L, and U
        [P(1, :), P(pivot_row, :)] = deal(P(pivot_row, :), P(1, :));
        [L(1, :), L(pivot_row, :)] = deal(L(pivot_row, :), L(1, :));
        [U(1, :), U(pivot_row, :)] = deal(U(pivot_row, :), U(1, :));
    end
    
    % Perform LU decomposition on the submatrix A(2:end, 2:end)
    temp = (U(2:end, 1) * U(1, 2:end)) / U(1, 1);
    [P_sub, L_sub, U_sub] = lu(U(2:end, 2:end) - temp);
    
    % Construct the final P, L, and U matrices
    P = P * [1, zeros(1,n-1); zeros(n-1,1), P_sub];
    L = [1, zeros(1,n-1); P_sub' * U(2:end,1) / U(1,1), L_sub];
    U = [U(1,1), U(1, 2:end); zeros(n-1,1), U_sub];
end