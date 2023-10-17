A = zeros(5);
A(1,2) = 1; A(2,3) = 1;
A(2,5) = 1; A(3,4) = 1;
A = 100*eye(5) + A + A';

L = rchol(A);

for i = 1:size(L, 1)
    fprintf('[');
    for j = 1:size(L, 2)
        fprintf('% .5f', L(i, j));
        if j < size(L, 2)
            fprintf(', ');
        end
    end
    fprintf(']\n');
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
