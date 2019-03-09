x = ptd_510050;
y = ptd_601668;

[ff, df] = initialize_exp_model(x, y)

function [ff, df] = initialize_exp_model(x, y)

    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    b = [2.5, -6, 2.5];
    [beta, R, ~, ~, ~, ~] = nlinfit(x, y, f, b);
    syms x;
    ff = beta(1) * exp(-beta(2) * x) + beta(3);
    df = diff(ff);

end