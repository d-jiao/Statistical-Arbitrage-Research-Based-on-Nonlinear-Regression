%clear, clc, warning off 
%fx=@(b,x)(b(1)+b(2)*x+b(3)*x.^2+b(4)*x.^3)./(1+b(5)*exp(-b(6)*x+b(7)*x.^2+b(8)*x.^3));
f = @(b, x) b(3)./(1 + (b(3) - b(1)) .* exp(b(2) .* (x)) ./ b(1));
figure(1), clf
x = ptd_601668;
y = ptd_510900;
plot(x, y, '.', 'markerfacecolor', 'b')
hold on
b = [0, -0.3, 1.1];
[beta, R, ~, CovB, MSE, ~] = nlinfit(x, y, f, b)
n = length(y);
SST = var(y) * (n - 1) %Total Sum of Squares
y1 = f(beta, x);
SSR = (y - y1)' * (y - y1) %Residual Sum of Squares
Rsquare = (SST - SSR) / SST
x1 = linspace(min(x), max(x), 350);
y1 = f(beta, x1);
plot(x1, y1, 'b-')
xlabel('601668.SH')
ylabel('510900.SH')
% axis tight
% legend('data1','fit1','data2','fit2','data3','fit3','location','best')