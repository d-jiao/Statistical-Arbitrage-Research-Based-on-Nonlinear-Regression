function QR = get_QR(x, y)
    [~, ~, res, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
    R = std(res);
    
    n = length(x) / 2;
    coeff = [];
    for i = 1 : n
        [beta, ~, ~, ~, ~] = regress(y(i + 1 : n + i), [ones(n, 1), x(i + 1 : n + i), x(i + 1 : n + i) .^ 2]);
        coeff = [coeff, beta];
    end
    coeff = coeff';
    Q = cov(diff(coeff));

    QR = [Q R];
end