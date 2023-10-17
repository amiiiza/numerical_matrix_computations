n = 100;

[Q, ~] = qr(rand(3));
A = Q * diag([10, 0.5, 0.1]) * Q';

[V, D] = eig(A);

eigenvalues = diag(D);

eigenvectors = V;

b = rand(3, 1);

alpha1 = zeros(1, n); 
alpha2 = zeros(1, n);
alpha3 = zeros(1, n);

z = zeros(3, n); 

z(:, 1) = b;

alpha1(1) = eigenvectors(:, 1)' * z(:, 1);
alpha2(1) = eigenvectors(:, 2)' * z(:, 1);
alpha3(1) = eigenvectors(:, 3)' * z(:, 1);

for i = 2:n 
    z(:, i) = A * z(:, i - 1);
    alpha1(i) = eigenvectors(:, 1)' * z(:, i);
    alpha2(i) = eigenvectors(:, 2)' * z(:, i);
    alpha3(i) = eigenvectors(:, 3)' * z(:, i);
end

disp('Eigenvalues:');
disp(eigenvalues);
disp('Alpha1:');
disp(alpha1);
disp('Alpha2:');
disp(alpha2);
disp('Alpha3:');
disp(alpha3);

iteration_steps = 1:n;

figure;
semilogy(iteration_steps, abs(alpha1), 'r', 'LineWidth', 1.5);
hold on;
semilogy(iteration_steps, abs(alpha2), 'g', 'LineWidth', 1.5);
semilogy(iteration_steps, abs(alpha3), 'b', 'LineWidth', 1.5);
title('Semilogarithmic Plot of |αij|');
xlabel('Iteration');
ylabel('|αij| (log scale)');
legend('|α1|', '|α2|', '|α3|');
grid on;