x = ptd_510050;
y = ptd_601668;
f = @(b, x) b(1) * exp(-b(2) * x) + b(3);

% l_m = @(x) b(1) + b(2) .* x + b(3) .* x .^ 2;

% f = @(b, x) b(3)./(1 + (b(3) - b(1)) .* exp(b(2) .* (x)) ./ b(1));
% x = ptd_601668;
% y = ptd_510900;

f1 = @(x) f(beta, x)

plot(x, y, '.', 'markerfacecolor', 'b')
hold on
plot(test_pthh(:,2), test_pthh(:,1), '.', 'markerfacecolor', 'r')
plot(test_pthh(:,2), f1(test_pthh(:,2)), '.', 'markerfacecolor', 'k')