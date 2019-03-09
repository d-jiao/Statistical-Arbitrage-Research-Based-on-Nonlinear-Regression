function [Q, R] = get_QR_exp(x, y)
    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    b = [2.5, -6, 2.5];
    [beta, res, ~, ~, ~, ~] = nlinfit(x, y, f, b);
    R = std(res);
    n = length(x) / 2;
    coeff = [];
    for i = 1 : n
        f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
        for j = 1 : 3
            [beta, ~, ~, ~, ~, ~] = nlinfit(x(i : i + n), y(i : i + n), f, beta);
            beta = lsqcurvefit(f, beta, x(i : i + n), y(i : i + n));
        end
        coeff = [coeff; beta];
    end
    Q = cov(diff(coeff));
end