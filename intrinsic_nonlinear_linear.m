f1 = @(x, eps) 1.7 * exp(0.7 .* x) .* eps
f2 = @(x, eps) 1.7 * exp(0.7 .* x) + eps
x = [0 : 0.2 : 10];
eps = randn(1, length(x));
y1 = f1(x, eps);
y2 = f2(x, eps);
plot(x, y1)
hold on
plot(x, y2)
legend('ģ��1','ģ��2��ģ��3', 'Location', 'south')