f = @(x, eps) 3 ./ (1 - 2 .* exp(-0.2 .* x) ./ 5) + eps
x = [0 : 0.2 : 50];
eps = randn(1, length(x))/50;
y = f(x, eps);
plot(x, y)