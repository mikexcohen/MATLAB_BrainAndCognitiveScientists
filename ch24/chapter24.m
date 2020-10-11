%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 24
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% simple 2D categorization

% use data from previous chapter
load ../ch23/times_090425blk10_ch115.mat

spikefeat = zeros(size(spikes,1),2);

for spikei=1:size(spikes,1)
    
    % find peak
    [~,peakidx] = max(spikes(spikei,:));
    
    % find minimum before peak
    [~,min1idx] = min(spikes(spikei,1:peakidx));
    
    % find minimum after peak
    [~,min2idx] = min(spikes(spikei,peakidx:end));
    min2idx = min2idx+peakidx-1;
    
    % get premin-peak difference
    spikefeat(spikei,1) = diff( spikes(spikei,[min1idx peakidx]) );
    
    % get min2min
    spikefeat(spikei,2) = min2idx-min1idx;
end

% plot one feature by the other
figure(1), clf
h = plot(spikefeat(:,1),spikefeat(:,2),'ko','markerfacecolor','k');
set(gca,'ylim',[2 12],'xlim',[-5 300])
set(h,'markeredgecolor',[.7 .7 .7])
xlabel('Spike amplitude')
ylabel('Spike width')

%% now with PCA

% subtract mean
spikes = bsxfun(@minus,spikes,mean(spikes,2));

% covariance matrix
spikecov = spikes'*spikes / (size(spikes,1)-1);

[eigvects,eigvals] = eig(spikecov);

figure(2), clf
subplot(221)
plot(mean(spikes,1))
set(gca,'xlim',[.75 12.25],'ylim',[-45 85])
xlabel('Time (a.u.)')
title('All spikes')

subplot(223)
plot(eigvects(:,end-1:end),'linew',2)
legend({'PC1';'PC2'})
set(gca,'xlim',[.75 12.25],'ylim',[-.7 .7])
title('Component time courses')
xlabel('Time (a.u.)')

subplot(222)
imagesc(spikecov)
set(gca,'clim',[-3000 3000])
xlabel('Time point'), ylabel('Time point')
title('Covariance matrix')

%% now with ICA

r = rank(spikecov);

weights = jader(spikes',r);
icas = weights*spikes';

figure(3), clf
plot(icas(1,:),icas(2,:),'k.','markersize',5)
xlabel('IC #1'), ylabel('IC #2')
set(gca,'xlim',[-7 0],'ylim',[-6 6])

% for Octave, kmeans is in the statistics package
ICclustidx = kmeans(icas(1:2,:)',2);

figure(4), clf
plot(icas(1,ICclustidx==1),icas(2,ICclustidx==1),'r.','markersize',3), hold on
plot(icas(1,ICclustidx==2),icas(2,ICclustidx==2),'b.','markersize',3)
xlabel('IC #1'), ylabel('IC #2')
legend({'cluster 1';'cluster 2'})

%% project each spike onto first and second components

comp1 = spikes*eigvects(:,end);
comp2 = spikes*eigvects(:,end-1);


figure(5)
plot(comp1,comp2,'k.','markersize',5)
xlabel('PC #1'), ylabel('PC #2')
% The axis limits were adjusted for the book, although it cuts off a few
% spikes. same for the next figure, and for the ICA-based analyses.
set(gca,'xlim',[-280 -0],'ylim',[-120 120])
%

% separate into two groups (you'll learn about kmeans in chapter 31).
PCclustidx = kmeans([comp1 comp2],2);

figure(6), clf
h1 = plot(comp1(PCclustidx==1),comp2(PCclustidx==1),'r^','markersize',4); hold on
h2 = plot(comp1(PCclustidx==2),comp2(PCclustidx==2),'bo','markersize',4);
xlabel('PC #1'), ylabel('PC #2')
set(gca,'xlim',[-280 -0],'ylim',[-120 120])
legend({'cluster 1';'cluster 2'})

set(h1,'color','k')
set(h2,'color',ones(1,3)*.5)

%% end.
