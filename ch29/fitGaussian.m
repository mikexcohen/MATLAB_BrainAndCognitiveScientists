function [sse,outparams] = fitGaussian(params,x,y)

%% setup Gaussian
peak = params(1); % amplitude of gaussian at peak
fwhm = params(2); % width (aka standard deviation)
cent = params(3); % x-axis point of peak lcation


% predicted value
predY = peak * exp( (-(x-cent).^2) ./ (2*fwhm^2) );

% predictors and sse (minimization objective)
outparams = [ peak fwhm cent ];
sse = sum( (predY-y).^2 ) / sum( y.^2 );

% optional plotting
plot(x,y,'ro',x,predY,'k-'); drawnow; pause(.01)

%%
