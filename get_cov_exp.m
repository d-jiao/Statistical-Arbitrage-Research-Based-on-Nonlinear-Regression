x = ptd_510050;
y = ptd_601668;
n = length(x) / 2;

f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
b = [2.5, -6, 2.5];
[beta, res, ~, ~, ~, ~] = nlinfit(x, y, f, b);
R = std(res);

coeff = [];

for i = 1 : n
    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    for j = 1 : 5
        [beta, ~, ~, ~, ~, ~] = nlinfit(x(i : i + n), y(i : i + n), f, beta);
        beta = lsqcurvefit(f, beta, x(i : i + n), y(i : i + n));
    end
    coeff = [coeff; beta];
end

Q = cov(diff(coeff));