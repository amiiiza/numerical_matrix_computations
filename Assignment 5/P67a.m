L1 = [11, 20, 100];
L2 = [10, 10, 10];
B = rand(2);
[Q, R] = qr(B);
b = [1;1];

for i = 1:length(L1)
    A = Q' * [L1(i) 0; 0 L2(i)] * Q;

    [x, y] = meshgrid(-2:0.1:2, -2:0.1:2);
    J_values = zeros(size(x));

    for row = 1:size(x, 1)
        for col = 1:size(x, 2)
            J_values(row, col) = 0.5 * ([x(row, col); y(row, col)]' * A * [x(row, col); y(row, col)]) - b'*[x(row, col); y(row, col)];
        end
    end

    figure;
    contour(x, y, J_values, 50);
    title(['Contour Plot for Condition Number = ', num2str(L1(i)/L2(i))]);
    xlabel('x');
    ylabel('y');
    colorbar;
end
