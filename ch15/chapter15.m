%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 15
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% theoretical and empirical distributions

figure(1), clf

obsval = 1.8;

% plot "theoretical" distribution
subplot(121)
x=-3:.001:3;
plot(x,exp(-x.^2))
hold on
plot([obsval obsval],get(gca,'ylim'),'m')
set(gca,'xlim',[-4 4])
xlabel('Statistic values')
ylabel('Probability')


% plot empirical distribution
% NOTE: You must first create Figure 15.3 for the next few lines to work.
subplot(122)
hist(permdiffs,50) 
hold on
plot([obsval obsval],get(gca,'ylim'),'m')
set(gca,'xlim',[-4 4])
xlabel('Statistic values')
ylabel('Count')

%% shifted and slightly skewed distribution

% create Gaussian distribution
r = randn(10000,1);

% get histogram values for later comparison
[y1,x1] = hist(r,100);


% skew the distribution
r(r>0) = log(1+r(r>0));
% and make it positive
r = -r-min(r);

% and get its histogram values
[y2,x2] = hist(r,100);



figure(2), clf
plot(x1,y1), hold on
plot(x2,y2,'r')
legend({'Gaussian';'positive skewed'})


%% simulated firing rate distributions

% number of trials
N = 100;

% male pictures
r = randn(N,1);
r(r>0) = log(1+r(r>0));
fr_males = 26-r*10;

% get histogram values for later comparison
[y1,x1] = hist(fr_males,20);


% female pictures
r = randn(N,1);
r(r>0) = log(1+r(r>0));
fr_females = 30-r*10;

% get histogram values for later comparison
[y2,x2] = hist(fr_females,20);


figure(3), clf
plot(x1,y1), hold on
plot(x2,y2,'r')
legend({'Males';'Females'})

%% mix trials together

% concatenate trials
allfr = cat(1,fr_males,fr_females);

% the following line does the same thing as the previous
% allfr = [fr_males fr_females];

% and mix them up (but in what order??)
allfr = allfr(randperm(N*2));


% a better approach
allfr = cat(1,fr_males,fr_females);

% repeat original and new order variables for accurate tracking
[conds,neworder] = deal( randperm(N*2) );

% re-sort data
allfr = allfr(neworder);

% and relable the conditions
conds(neworder<N+1) = 1;
conds(conds>1) = 0;

%% generate one null hypothesis scenario

% random permutation
fakeconds = randperm(N*2);

% shuffled condition labels
fakeconds(fakeconds<N+1) = 1;
fakeconds(fakeconds>1) = 0;


% these two means should be different.
[mean(allfr(conds==1)) mean(allfr(conds==0))]

% should these two be different?
[mean(allfr(fakeconds==1))  mean(allfr(fakeconds==0)) ]

%% and now a distribution of null hypothesis values

nPerms = 1000;
permdiffs = zeros(nPerms,1);

for permi=1:nPerms
    fconds = randperm(N*2);
    fconds(fconds<N+1) = 1;
    fconds(fconds>1) = 0;
    permdiffs(permi) = mean(allfr(fconds==0))-mean(allfr(fconds==1));
end


% plot the distribution of H0 values
figure(4), clf
hist(permdiffs,50)
hold on

% and plot the observed value on top
obsval = mean(allfr(conds==0))-mean(allfr(conds==1));
plot([obsval obsval],get(gca,'ylim'),'m','linew',10)
xlabel('Value'), ylabel('Count')

%% two methods of evaluating statistical significance

% Z-value
zVal = ( obsval-mean(permdiffs) ) / std(permdiffs);
p = 1-normcdf(abs(zVal));

% p-value count
pCount = sum(permdiffs>obsval)/nPerms;

%% now with real data

load tfdata

% p-value threshold
pval = .05;

% the empirically observed time-frequency power difference
realdif = squeeze( mean(allpow(:,:,101:200),3) - mean(allpow(:,:,1:100),3) );

% setup permutation testing
nPerms  = 1000; % number of iterations
permdif = zeros([ nPerms size(realdif) ]);

% loop over permutations
for permi=1:nPerms
    
    % generate a fake trial order
    fakeord = randperm(size(allpow,3));
    
    % compare this to how the variable 'realdif' was created
    fakedif = squeeze( mean(allpow(:,:,fakeord(1:100)),3) - mean(allpow(:,:,fakeord(101:200)),3) );
    
    % does this need to be a separate line from the previous?
    permdif(permi,:,:) = fakedif;
end

% compute z difference score
permmean = squeeze(mean(permdif,1));
permstd  = squeeze(std(permdif,[],1));
zmap = ( realdif-permmean ) ./ permstd;


% plotting...
figure(5), clf
subplot(221)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Power condition difference')


subplot(222)
contourf(timevec,frex,zmap,40,'linecolor','none')
colorval = max(abs(zmap(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Power condition Z-map')


subplot(223)

% Threshold the map by turning all subtreshold pixels to 0.
zthresh = zmap;
zthresh(abs(zthresh)<norminv(1-pval)) = 0;

contourf(timevec,frex,zthresh,40,'linecolor','none')
colorval = max(abs(zthresh(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Thresholded differences')


subplot(224)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Significant regions outlined')

% Draw significance contours on top of the colorful map.
hold on
contour(timevec,frex,logical(zthresh),1,'k')

%% extreme-value-based correction

% p-value threshold
pval = .05;

% initialize extreme-value distribution
exvals = zeros(nPerms,2);

for permi=1:nPerms
    
    % again with the shuffled trial order
    fakeord = randperm(size(allpow,3));    
    fakedif = squeeze( mean(allpow(:,:,fakeord(1:100)),3) - mean(allpow(:,:,fakeord(101:200)),3) );

    % from this map, take the most extreme positive and negative values
    exvals(permi,1) = min(fakedif(:));
    exvals(permi,2) = max(fakedif(:));
end

% find the thresholds
lowerThresh = prctile(exvals(:,1),100*pval);
upperThresh = prctile(exvals(:,2),100-100*pval);


% let's have a look at the distribution of extreme values
figure(6), clf
histogram(exvals(:),200,'DisplayStyle','stairs')
xlabel('Value'), ylabel('Count')

% plot the distribution of real values on top
hold on
h=histogram(realdif(:),200);

% and draw the thresholds
plot(ones(1,2)*lowerThresh,get(gca,'ylim'),'k--')
plot(ones(1,2)*upperThresh,get(gca,'ylim'),'k--')

% finally, change some of the histogram properties
set(h,'edgecolor','none','facecolor','r')

%%

% threshold the real map
threshmap = realdif;
threshmap(threshmap>lowerThresh & threshmap<upperThresh) = 0;



figure(7), clf
subplot(221)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])

subplot(222)
contourf(timevec,frex,threshmap,40,'linecolor','none')
colorval = max(abs(threshmap(:)))*.8;
set(gca,'clim',[-colorval colorval])

subplot(223)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])

hold on
contour(timevec,frex,logical(threshmap),1,'k')

%% end
