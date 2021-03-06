function [fitresult, gof] = createFit(tnew, quant_new)
%CREATEFIT(TNEW,QUANT_NEW)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tnew
%      Y Output: quant_new
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 25-Jul-2018 14:16:58


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( tnew, quant_new );

% Set up fittype and options.
ft = fittype( 'a*tanh((x-C)/tau)+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.45 13 10 1.95];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'quant new vs. tnew', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel tnew
ylabel quant_new
grid on
hold on 

