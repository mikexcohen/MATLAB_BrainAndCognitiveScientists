%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 2
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% Figure 2.01
% using polynomials to demonstrate fitting vs. overfitting
% (you'll learn about what this code means in chapter 28)

% generate data
N = 10; % number of datapoints
data1 = linspace(1,10,N) + 2*randn(1,N); % linearly increasing with noise
data2 = linspace(1,10,N) + 2*randn(1,N);

% polynomial that over-fits the data (N+1 coefficients!)
polycoefsOver = polyfit(1:N,data1,N);
% estimated data based on the coefficients
yHat = polyval(polycoefsOver,1:N);

% plot
figure(1), clf

% plot the data and over-fit model
subplot(221)
plot(yHat,'r','linew',4), hold on
plot(data1,'ko','linew',3,'markerface','y','markersize',15)
set(gca,'xlim',[0 11],'ylim',[-5 20])
title('Overfitting')

% plot the model on new data with the same underlying 
% structure but different noise.
subplot(223)
plot(yHat,'r','linew',4), hold on
plot(data2,'ko','linew',3,'markerface','y','markersize',15)
set(gca,'xlim',[0 11],'ylim',[-5 20])



% now a polynomial fit with 1 coefficient
polycoefsGood = polyfit(1:N,data1,1);
% estimated data based on the coefficients
yHat = polyval(polycoefsGood,1:N);

% plot data and model fit
subplot(222)
plot(yHat,'r','linew',4), hold on
plot(data1,'ko','linew',3,'markerface','y','markersize',15)
set(gca,'xlim',[0 11],'ylim',[-5 20])
title('Fitting')

% plot new data and same model
subplot(224)
plot(yHat,'r','linew',4), hold on
plot(data2,'ko','linew',3,'markerface','y','markersize',15)
set(gca,'xlim',[0 11],'ylim',[-5 20])

%% end.
