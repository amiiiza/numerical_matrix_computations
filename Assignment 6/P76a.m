n = 20;
[Q, R] = qr(rand(n));

i_values = 5:15;
L_values = [1, 10, 100, 1000];
b = rand(n, 1)*10;

for L_index = 1:length(L_values)
    L = L_values(L_index);
    A = Q * diag(linspace(1, L, n)) * Q';
    x = A \ b;
    errors = zeros(size(i_values));
    
    for i_index = 1:length(i_values)
        i = i_values(i_index);
        [Q_krylov, R_krylov] = my_arnoldi(A, b, i);
        M = Q_krylov' * A * Q_krylov;
        L_krylov = Q_krylov' * b;
        q_i = linsolve(M, L_krylov);
        x_i = Q_krylov * q_i;
        error = sqrt((x_i - x)'* A *(x_i-x));
        errors(i_index) = error;
    end

    figure;
    plot(i_values, errors, '-o', 'LineWidth', 1.5);
    xlabel('i');
    ylabel('Error (||x_i - x||_A');
    title(['Error Dependence on i (L = ', num2str(L), ')']);
end

function [Q, R] = my_arnoldi(A, b, N)
    Q = [];
    q = b;
    for i = 1:N
        for k = 1:size(Q, 2)
            R(k, i) = q' * Q(:, k);
            q = q - R(k, i) * Q(:, k);
        end
        R(i, i) = norm(q);
        Q(:, i) = q / R(i, i);
        q = A * Q(:, i);
    end
end