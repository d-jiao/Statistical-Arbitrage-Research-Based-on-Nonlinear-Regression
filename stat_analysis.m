%qqplot(R);

figure
subplot(2,1,1)
autocorr(R)
subplot(2,1,2)
parcorr(R)
% acf = autocorr(R,15);
% pacf = parcorr(R,15);
% [h, p, Qstat, crit] = lbqtest(R, 'Lags', [5, 10, 15]);

% figure
% subplot(2,1,1)
% autocorr(R.^2)
% subplot(2,1,2)
% parcorr(R.^2)

% [h, p, Qstat, crit] = lbqtest(R .^ 2, 'Lags', [5, 10, 15])

% LM test
% e = R - mean(R);
% [h,p,fStat,crit] = archtest(e,'Lags',2)