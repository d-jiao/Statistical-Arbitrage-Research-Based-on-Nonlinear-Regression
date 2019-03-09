f = @(x, eps) 5 ./ (1 + 2 .* exp(-0.2 .* x) ./ 3) + eps
x = [0 : 0.2 : 50];
eps = randn(1, length(x))/50;
y = f(x, eps);
plot(x, y)