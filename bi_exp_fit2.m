%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 18-May-2017 17:21:29 自动生成


%% Fit: 'untitled fit 1'.
x = ptd_600028;
y = ptd_510900;

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-2.55772375240779e-08 2.97271035847391 0.105698345461799 0.469268295896909];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

plot(xData, yData, '.', 'markerfacecolor', 'b')
hold on
x1 = linspace(min(x), max(x), 350);
y1 = fitresult(x1);
plot(x1, y1, 'b-')
xlabel('600028.SH')
ylabel('510900.SH')
% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, 'b' );
% %legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
%grid on

n = length(yData);
y1 = fitresult(x);
SST = var(yData) * (n - 1) %Total Sum of Squares
SSR = (yData - y1)' * (yData - y1) %Residual Sum of Squares
Rsquare = (SST - SSR) / SST


