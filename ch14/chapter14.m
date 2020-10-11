%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 14
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% d'

% step 1
hitP = 22/30;
faP  =  3/30;

% step 2
hitZ = norminv(hitP);
faZ  = norminv(faP);

% step 3
dPrime = hitZ-faZ;

%% aside on norminv

x = 0:.001:1;

p1 = norminv(x); % with stats toolbox
p2 = -sqrt(2).*erfcinv(2*x); % without stats toolbox, also works in Octave

figure(1), clf
plot(p1,x,'b','linew',8), hold on
plot(p2,x,'ro')
xlabel('Z value'), ylabel('proportion')

%% 2D d' space

% convert counts to probabilities
x  = .01:.01:.99;

dp = bsxfun(@minus,norminv(x)',norminv(x));

% show the 2D d' space
figure(2), clf
contourf(x,x,dp,80,'linecolor','none')
xlabel('False alarm rate')
ylabel('Hit rate')
hold on % we'll plot lines on top
axis square

% colors for isosensitivity curves
colorz = 'rbmk';

% the d' values
dp2plot = [1 1.5 2 2.5];
tol = .01;

for dpi=1:length(dp2plot)
    
    % find points
    idx = find(dp>dp2plot(dpi)-tol & dp<dp2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(dp),idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',4,'markersize',9)
end

%% response bias

% step 1
hitP = 22/30;
faP  =  3/30;

% step 2
hitZ = norminv(hitP);
faZ  = norminv(faP);

% step 3
dPrime = hitZ-faZ;
respBias = -(hitZ+faZ)/2;

%% 2D bias space

% convert counts to probabilities
rb = -bsxfun(@plus,norminv(x)',norminv(x))/2;


figure(3), clf
contourf(x,x,rb,80,'linecolor','none')
xlabel('False alarm rate')
ylabel('Hit rate')
hold on
axis square


colorz = 'rbmk';

rb2plot = [.3 .5 .9 1.5];
tol = .01;

for dpi=1:length(rb2plot)
    
    % find points
    idx = find(rb>rb2plot(dpi)-tol & rb<rb2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(rb),idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',4,'markersize',9)
end

%% discretization

% define parameters
ntrials = 100;
nbins   = 7;

% create random data
d = [500+100*randn(ntrials,1) rand(ntrials,1)>.3];
d = sortrows(d,1);

% discretization based on n bins
binidx = ceil(linspace(0,nbins-1,length(d)));
discdata = zeros(2,nbins);

for i=1:nbins
    discdata(1,i) = mean(d(binidx==i,1));
    discdata(2,i) = mean(d(binidx==i,2));
end

% and plot.
figure(4), clf
plot(discdata(1,:),discdata(2,:),'o-','markerfacecolor','r','markersize',10)

%% discretization using tiedrank

% more random data
d = [500+100*randn(ntrials,1) rand(ntrials,1)>.3];

% procedure in separate steps
temp  = tiedrank(d(:,1))/ntrials; % tiedrank is not in default Octave, but you can find it online
temp  = nbins*temp;
drank = ceil( temp );

% can also be done in one line
drank = ceil( nbins*tiedrank(d(:,1))/ntrials );

%% CAFs with real data

load behavioralDataRK.mat

% number of discretizations
nbins = 12;

% and now discretize
drank = ceil( nbins*tiedrank(beh(:,2))/length(beh) );

% initialize and compute
caf = zeros(nbins,2);
for i=1:nbins
    caf(i,1) = mean(beh(drank==i,2));
    caf(i,2) = mean(beh(drank==i,1));
end

% and plot
figure(5), clf
plot(caf(:,1),caf(:,2),'o-','markerfacecolor','r','markersize',10)
xlabel('Average RT')
ylabel('Average accuracy')

%% end
