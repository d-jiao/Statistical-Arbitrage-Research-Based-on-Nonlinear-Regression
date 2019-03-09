function max_drawdown = max_drawdown(dta)
    max_drawdown = 0;
    n = length(dta);
    for i = 1 : n - 2
        for j = i + 1 : n
            max_drawdown = min(max_drawdown, dta(j) - dta(i));
        end
    end
end