n_values = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
expression_values = [];
gamma_values = [];
decimal_places = 2;

for i = 1:14
    k = n_values(i);
    
    n = 2^k;

    a21 = zeros(n-1, 1);

    for j = 1:n-1
        a21(j, 1) = sqrt(pi + j);
    end
    
    a21_sq = (n-1) * pi + (n * (n - 1)) / 2;

    a11 = a21_sq + 1;

    I = speye(n - 1);

    A = sparse([a11, a21'; a21, I]);

    e1 = zeros(1, n);
    e1(1) = 1;
  
    x = e1 / A;

    actual = [1, -a21'];

    error = abs(x - actual);
    
    error_ratio = error ./ actual;

    expression_values = [expression_values,max(error_ratio(:))];

end

figure;
semilogy(n_values, expression_values, 'b-o', 'LineWidth', 2, 'DisplayName', 'relative error');
hold on;
ylabel('Value (log scale)');
title('Relative error of actual value and matlab value');
legend('Location', 'Best');
grid on;


