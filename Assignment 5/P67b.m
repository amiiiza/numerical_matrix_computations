max_iterations = 4000;
threshold = eps;
L1 = [110, 120, 100];
L2 = [1, 1, 1];

for i = 1:length(L1)
     B = rand(2);
     [Q,R] = qr(B);
     A = Q' * [L1(i) 0 ; 0 L2(i)] * Q;
     x_answer = [rand(1); rand(1)];
     b = A * x_answer;

    x_current = [0; 0];
    error_0 = (x_answer - x_current)' * A * (x_answer - x_current);
    E = [];
    predicted_error_rate = [];
    for iter = 1:max_iterations
        p = b - A * x_current;
        r = b - A * x_current;
        alpha = (p' * r) / (p' * A * p);
        
        x_current = x_current + alpha * p;

        error = x_answer - x_current;
        disp((1 - (cond(A)^(-2))));
        E = [E, error' * A * error];
        predicted_error_rate = [predicted_error_rate,(1 - (cond(A)^(-2))^iter)];
        if (error' * A * error) < threshold
           disp(x_answer);
           disp(x_current);
           break
        end
    end

    figure;
    semilogy(E);
    title(['Error in A-weighted ', i]);
    xlabel('Iteration');
    ylabel('Error (log scale)');
    hold on;
    semilogy(predicted_error_rate);
    legend('Actual Error', 'Predicted Error Rate');
    hold off;
end
