f = @(b, x) b(1) * exp(-b(2) .* x) + b(3);
%f = @(b, x) b(1) * exp(-b(2) * (x - b(3)));
figure(1), clf
x = ptd_000839;
y = ptd_601668;
plot(x, y, '.', 'markerfacecolor', 'b')
hold on
b = [70, 0.5, 5];
for l=1:5
    b=lsqcurvefit(f,b,x,y);
    [b, R, ~ , Cov, MSE, ~]=nlinfit(x,y,f,b);
end
% [beta, R, ~, CovB, MSE, ~] = nlinfit(x, y, f, b);
n = length(y);
SST = var(y) * (n - 1); %Total Sum of Squares
y1 = f(b, x);
SSR = (y - y1)' * (y - y1); %Residual Sum of Squares
Rsquare = (SST - SSR) / SST;
x1 = linspace(min(x), max(x), 350);
y1 = f(b, x1);
plot(x1, y1, 'b-')
xlabel('000839.SH')
ylabel('601668.SH')
%axis tight