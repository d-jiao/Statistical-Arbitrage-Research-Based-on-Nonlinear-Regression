f = @(x, eps) 35 * exp(-0.7 .* x) + eps
x = [0 : 0.05 : 10];
eps = randn(1, length(x));
y = f(x, eps);
plot(x, y)