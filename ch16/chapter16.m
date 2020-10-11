%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 16
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% generate bivariate correlated data

r = .6;
n = 100;
x = randn(n,2);
x(:,2) = x(:,1)*r + x(:,2)*sqrt(1-r^2);

figure(1), clf
plot(x(:,1),x(:,2),'o')
obsR = corr(x);
title([ 'R=' num2str(obsR(2)) ])

%% multivariate covariance with real data

load ../sampleEEGdata

% extract channels-by-time matrix from trial 1
data = squeeze(EEG.data(:,:,1));

% is the a valid covariance?
cov2 = data*data'/(EEG.pnts-1);

% this one's better (why?)
data = bsxfun(@minus,data,mean(data,2));
cov2 = data*data'/(EEG.pnts-1);


% show the covariance matrix
figure(2), clf
imagesc(cov2)
set(gca,'clim',[-100 200])


% now show the variance ("co"variance along the diagonal)
figure(3), clf
topoplotIndie(diag(cov2),EEG.chanlocs);
set(gca,'clim',[0 200])
colormap hot
colorbar
title('Variances')


% and some seeded covariances. The numbers in cov2
% are particular electrodes.
figure(4), clf
subplot(131)
topoplotIndie(cov2(47,:),EEG.chanlocs);

subplot(132)
topoplotIndie(cov2(23,:),EEG.chanlocs);

subplot(133)
topoplotIndie(cov2(10,:),EEG.chanlocs);

%% correlation in one variable

r = .6;
n = 100;
x = randn(n,2);
x(:,2) = x(:,1)*r + x(:,2)*sqrt(1-r^2);

% stretch to make covariance different from correlation
x(:,1) = x(:,1)*100 + 2;
x(:,2) = x(:,2)*20 + 10;

cx  = x(:,1)'*x(:,2) / (n-1);
vx1 = x(:,1)'*x(:,1) / (n-1);
vx2 = x(:,2)'*x(:,2) / (n-1);

cov1 = cov(x)

% these two variables should be identical!
cor1 = cx / sqrt(vx1*vx2)
cor2 = corr(x)

%% correlation

% extract and prepare the data
data = squeeze(EEG.data(:,:,1));
data = bsxfun(@minus,data,mean(data,2));
data = bsxfun(@rdivide,data,std(data,[],2));

% compute correlation as covariance of normalized data
cor2 = data*data'/(EEG.pnts-1);

% image it
figure(5), clf
imagesc(cor2)
get(gca,'clim') % let's see the color range

%% or...

% data again
data = squeeze(EEG.data(:,:,1));
data = bsxfun(@minus,data,mean(data,2));

% this is the unscaled covariance
cor2 = data*data'/(EEG.pnts-1);


% compute the variance matrix
stdMat = sqrt( diag(cor2)*diag(cor2)' );

% and then scale the (previously unscaled) covariance
cor2 = cor2 ./ stdMat;

imagesc(cor2)
get(gca,'clim') % how does this color range compare to that of the previous cell?

%% Anscobe's quartet to illustrate Pearson vs. Spearman

anscombe = [
  % series 1    series 2    series 3     series 4
    10  8.04    10  9.14    10  7.46     8  6.58;
     8  6.95     8  8.14     8  6.77     8  5.76;
    13  7.58    13  8.76    13 12.74     8  7.71;
     9  8.81     9  8.77     9  7.11     8  8.84;
    11  8.33    11  9.26    11  7.81     8  8.47;
    14  9.96    14  8.10    14  8.84     8  7.04;
     6  7.24     6  6.13     6  6.08     8  5.25;
     4  4.26     4  3.10     4  5.39     8  5.56;
    12 10.84    12  9.13    12  8.15     8  7.91;
     7  4.82     7  7.26     7  6.42     8  6.89;
     5  5.68     5  4.74     5  5.73    19 12.50;
    ];


%% end
