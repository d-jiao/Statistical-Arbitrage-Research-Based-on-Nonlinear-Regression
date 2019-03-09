x = ptd_510050;
y = ptd_601668;
[b, bint, r, rint, stats] = regress(y, [ones(size(x)), x, x .^ 2]);
l_m = @(x) b(1) + b(2) .* x + b(3) .* x .^ 2;
y1 = l_m(x);
n = length(y);
SST = var(y) * (n - 1) %Total Sum of Squares
y1 = l_m(x);
SSR = (y - y1)' * (y - y1) %Residual Sum of Squares
Rsquare = (SST - SSR) / SST
x1 = linspace(min(x), max(x), 350);
y1 = l_m(x1);
plot(x, y, '.', 'markerfacecolor', 'b')
hold on
plot(x1, y1, 'b-')
xlabel('510050.SH')
ylabel('601668.SH')