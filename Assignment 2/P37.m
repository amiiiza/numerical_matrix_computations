A2 = [
    20 0 1 1 1 1 0;
    0 20 1 1 0 0 1;
    1 1 20 0 0 0 0;
    1 1 0 20 0 0 0;
    1 0 0 0 20 0 0;
    1 0 0 0 0 20 0;
    0 1 0 0 0 0 20;
];

P = my_md(A2)';
fprintf('a)\n');
disp(P);

L = rchol(A2);
LP = rchol(P'* A2 * P);
numNonZeros_L = nnz(L);
numNonZeros_PTA2P = nnz(LP);

fprintf('b) Number of non-zeros in Cholesky factor of A2: %d\n', numNonZeros_L);
fprintf('Number of non-zeros in Cholesky factor of PT A2P: %d\n', numNonZeros_PTA2P);

fprintf('\nc)\n');
P_amd = amd(A2)';

disp(P_amd);

L = rchol(A2);
LP = rchol(P_amd'* A2 * P_amd);
numNonZeros_L = nnz(L);
numNonZeros_P_amdTA2P_amd = nnz(LP);

fprintf('Number of non-zeros in Cholesky factor of A2: %d\n', numNonZeros_A2);
fprintf('Number of non-zeros in Cholesky factor of P_amdT A2P_amd: %d\n', numNonZeros_P_amdTA2P_amd);


function p = my_md(A)
    n = size(A, 1);
    p = 1:n;
    
    for i = 1:(n - 1)
        % Try all remaining entries as entry i
        nnzLi = zeros(1, n);
        
        for j = (i + 1):n
            tmp = p;
            tmp(i) = p(j);
            tmp(j) = p(i);
            nnzLi(j) = length(unique(my_reach(A(tmp, tmp), i, 1:(i - 1))));
        end
        
        % Choose permutation minimizing nnz in column i.
        [~, I] = min(nnzLi((i + 1):n));
        I = I(1) + i;
        pi = p(i);
        p(i) = p(I(1));
        p(I) = pi;
    end
end

function [R, visited] = my_reach(A, v, S, R, visited)
    if nargin == 3
        R = [];
        visited(1:size(A, 2)) = false;
    end
    visited(v) = true;
    edges = find(abs(A(:, v)) > 0);
    if isempty(S)
        R = setdiff(edges, v);
        return;
    end
    for w = edges(:)'
        if ~visited(w)
            if ~ismember(S, w)
                R = [R w];
            else
                [R, visited] = my_reach(A, w, S, R, visited);
            end
        end
    end
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

