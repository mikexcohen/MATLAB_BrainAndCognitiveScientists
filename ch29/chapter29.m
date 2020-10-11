%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 29
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% fitting a simple function

% comment one of the lines below

funch = @(t) 1*(t).^2 + 0*t + 0; % try changing the parameters (book sets several implicitly to 0)

[xval,funcval] = fminsearch(funch,-2);

figure(1), clf

x = -4:.001:4;
plot(x,funch(x),'k','markersize',5), hold on
plot(xval,funcval,'ro','markerface','r'), hold off
axis square
xlabel('X'), ylabel('Y')

%% fitting a piece-wise linear model

% generate a triangle distribution
a = .2;
c = .6;
b = .9;
x = rand(1,10000);

y(x<c) = a + sqrt( x(x<c).*(b-a).*(c-a) );
y(x>c) = b - sqrt( (1-x(x>c)).*(b-a).*(b-c) );

% convert x and y into a distribution
[y,x] = hist(y,100);

% find a decent starting location
[~,initB] = min(abs(x-.5));

% create a function handle
funch = @(initB) fit2segLinear(initB,x,y);

% and fit the model!
figure(2), clf
[optBreakPoint,sse,exitflag,fmininfo] = fminsearch(funch,initB);
% octave note: fminsearch only provides the first two outputs

%% now with gaussian

% parameters of the Gaussian
peak = 4;
fwhm = 1;
cent = 3;
nois = .5;

x = -10:.1:10;

gaus = peak*exp( -(x-cent).^2 / (2*fwhm^2) );
gaus = gaus + nois*randn(size(gaus));


% initialize:  peak  fwhm  center
initParms = [   2     2      -2   ];
funch = @(initParms) fitGaussian(initParms,x,gaus);

figure(3), clf
[outparams,sse,exitflag,fmininfo] = fminsearch(funch,initParms);

%% show getting stuck in a local minimum using an inverse sinc function

% define sinc function and handle
funch = @(x) -sin(x)./x;

% minimize function
[xval,funcval] = fminsearch(funch,0);

figure(5), clf
x = -50:.01:100;
plot(x,funch(x))
hold on
plot(xval,funcval,'ro','markersize',10,'markerfacecolor','k')

%% hist function

x = randn(1000,1);

figure(6), clf
hist(x,40)

% get outputs and it doesn't plot
[yy,xx] = hist(x,40);

hold on
plot(xx,yy,'r','linew',3)

%% histogram

x = randn(1000,1);

figure(7), clf
histogram(x,40) % this function doesn't exist (yet) in octave

% get outputs and it doesn't plot
hdata = histogram(x,40);

hold on
xvals = ( hdata.BinEdges(1:end-1) + hdata.BinEdges(2:end) )/2;
plot(xvals,hdata.Values,'r','linew',3)

%% assigning data points to bins

xidx = zeros(size(x));

for bini=1:hdata.NumBins
    ix = hdata.Data > hdata.BinEdges(bini) & hdata.Data < hdata.BinEdges(bini+1);
    xidx(ix) = bini;
end

%% end.
