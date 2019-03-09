function sharpe_ratio = sharpe_ratio(dta)
    rtd = diff(log(dta));
    sharpe_ratio = mean(rtd) / std(rtd);
end