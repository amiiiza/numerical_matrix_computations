n = 5;       

a21 = zeros(n-1, 1);

for i = 1:n-1
    a21(i, 1) = sqrt(pi + i);
end

a21_sq = (n-1) * pi + (n * (n - 1)) / 2;

a11 = a21_sq + 1;
I = speye(n - 1);

A = sparse([a11, a21'; a21, I]);

disp(A)