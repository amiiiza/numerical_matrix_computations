n_values = [10, 20, 40, 80, 160, 320];
expression_values = [];
gamma_values = [];
decimal_places = 2;

for i = 1:6
    n = n_values(i);
    
    a21 = zeros(n-1, 1);

    for j = 1:n-1
        a21(j, 1) = sqrt(pi + j);
    end
    
    a21_sq = (n-1) * pi + (n * (n - 1)) / 2;

    a11 = a21_sq + 1;

    I = speye(n - 1);

    A = sparse([a11, a21'; a21, I]);
  
    Lb = chol(A, 'lower');

    error = abs(A - Lb * Lb');
    
    product = abs(Lb) * abs(Lb');
    
    error_ratio = error ./ product;

    expression_values = [expression_values,max(error_ratio(:))];

    gamma_values = [gamma_values,((n+1) * eps)/(1 - (n+2) * eps)];

end

figure;
semilogy(n_values, expression_values, 'b-o', 'LineWidth', 2, 'DisplayName', 'max\_ij |a_{ij} - (LLT)_{ij}| / (|L||LT|)_ij');
hold on;
semilogy(n_values, gamma_values, 'r--', 'LineWidth', 2, 'DisplayName', '\gamma_{n}');
xlabel('n');
ylabel('Value (log scale)');
title('Comparison of Expression and \gamma_{n}');
legend('Location', 'Best');
grid on;


