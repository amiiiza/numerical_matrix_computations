max_iterations = 6;
threshold = eps;
L1 = [11, 20, 100];
L2 = [100, 1000, 1000];

for i = 1:length(L1)
    B = rand(2);
    [Q,R] = qr(B);
    A = Q' * [L1(i) 0 ; 0 L2(i)] * Q;
    x_answer = [rand(1); rand(1)];
    b = A * x_answer;
    
    x_current = [0; 0];
    p = [];
    for iter = 1:max_iterations
        beta = [];
        r = b - A * x_current;

        for s = 1:iter-1
            t = [p(1,s);p(2,s)];
            beta = [beta, (r' * A * t)/(t' * A * t)];
            disp("Validating Beta's: ");
            disp(beta);
        end
        
        p_new = r;
        for s = 1:iter-1
            t = [p(1,s);p(2,s)];
            p_new = p_new - beta(s) * t;
        end

        p = [p, p_new];
        
        alpha = (p_new' * r) / (p_new' * A * p_new);
        
        x_current = x_current + alpha * p_new;

        if(norm(r) < eps)
            break;
        end
    end
end
