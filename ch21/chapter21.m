%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 21
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% basic RMS on empirical data

load eeg1

% compute RMS. pretty simple.
rmsx = sqrt( mean( eeg.^2 ,1) );

% plot the time-domain average response
figure(1), clf
subplot(211)
plot(EEGtime,mean(eeg,2))
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'xlim',EEGtime([1 end]))

% bar plot of the RMS per trial. Notice anything striking?
subplot(212)
bar(rmsx,'histc')
ylabel('RMS')
xlabel('Trials')
set(gca,'xlim',[0 size(eeg,2)+2])

%% RMS over space per time point

load ../sampleEEGdata

% compute RMS per trial or on trial-average
rmsx_trials = squeeze( sqrt( mean( EEG.data.^2 ,1) ) );
rmsx_averag = squeeze( sqrt( mean( mean(EEG.data,3).^2 ,1) ) );

% and plot
figure(2), clf
plot(EEG.times,mean(rmsx_trials,2),EEG.times,rmsx_averag)
legend({'All trials';'Trial average'})
xlabel('Time (ms)'), ylabel('\muV^2')

%% DFA on human task data

% load data
load behaviorData.mat

N = length(perfdat);

% setup parameters
nScales = 20;
ranges  = [1 400]*srate;
scales  = ceil(logspace(log10(ranges(1)),log10(ranges(2)),nScales));
rmses   = zeros(size(scales));

% plot the data
figure(3), clf
subplot(211)
plot(timevec/60,perfdat)
xlabel('Time (minutes)'), ylabel('Performance')
set(gca,'ytick',[-1 1],'yticklabel',{'error';'correct'},'ylim',[-1.1 1.1],'xlim',timevec([1 end])/60)

% integrate data
perfdat = cumsum(perfdat(:)-mean(perfdat));

% and show that for comparison
subplot(212)
plot(timevec/60,perfdat)
set(gca,'xlim',timevec([1 end])/60)
xlabel('Time (minutes)')


% qualitative demonstration of self-similar behavior
figure(4), clf
for i=1:4
    subplot(4,1,i)
    plot(timevec/60,detrend(perfdat))
end
%
subplot(411), set(gca,'xlim',[0 40])
subplot(412), set(gca,'xlim',[25 29])
subplot(413), set(gca,'xlim',[27 27.4])
subplot(414), set(gca,'xlim',[27.2 27.24])

%% an aside on detrend...

% 200 points, because why not?
n = 200;

% random noise plus a linear trend
d = randn(1,n) + linspace(-3,10,n);

figure(5), clf
plot(1:n,d, 1:n, detrend(d) )
legend({'Original';'Detrended'})

%% compute RMS over different time scales

%
for scalei = 1:length(scales)
    
    % epoch
    n = floor(N/scales(scalei)); % number of epochs
    epochs = reshape( perfdat(1:n*scales(scalei)) ,scales(scalei),n)';
    
    % detrend
    depochs = detrend(epochs')';
    % note the effect of detrending:
    %clf, subplot(211), plot(epochs'), subplot(212), plot(depochs')
    
    % RMS and average
    rmses(scalei) = mean( sqrt( mean(depochs.^2,1) ) );
end


% fit a linear model to quantify scaling exponent
A = [ ones(length(scales),1) log10(scales)' ];
dfa = (A'*A) \ (A'*log10(rmses)');

% plot the 'linear' fit (in log-log space)
figure(6), clf
plot(log10(scales/srate),log10(rmses),'ks','linew',3,'markerfacecolor','w','markersize',18)

hold on
plot(log10(scales/srate),dfa(1)+dfa(2)*log10(scales),'k--','linew',2)

%% now repeat the procedure using DMA

for scalei = 1:length(scales)
    
    % create kernel for this scale
    nConv  = N+scales(scalei)-1;
    kernel = fft(ones(scales(scalei),1)/scales(scalei),nConv);
    hfKrn  = floor(scales(scalei)/2)+1; % half of kernel
    
    % mean-smooth as convolution
    convres = ifft( fft(perfdat,nConv) .* kernel );
    convres = convres(hfKrn:end-hfKrn+1+mod(nConv,2));
    
    residX = perfdat - convres;
    
    
    n = floor(N/scales(scalei)); % number of epochs
    epochs = reshape(residX(1:n*scales(scalei)),scales(scalei),n)';
    rmses(scalei) = mean( sqrt( mean(epochs.^2,1) ) );
end

dma = (A'*A) \ (A'*log10(rmses)');


h1 = plot(log10(scales/srate),log10(rmses),'kd','linew',3,'markerfacecolor','w','markersize',18);

hold on
h2 = plot(log10(scales/srate),dma(1)+dma(2)*log10(scales),'k--','linew',2);

legend({'Data (DFA)';[ 'Fit (DFA m=' num2str(dfa(2)) ')' ];...
        'Data (DMA)';[ 'Fit (DMA m=' num2str(dma(2)) ')' ]})
xlabel('log10(Scales)')
ylabel('log10(RMS)')

set(h1,'color',ones(1,3)*.6)
set(h2,'color',ones(1,3)*.6)
set(gca,'xlim',[-.1 2.8],'ylim',[1.2 4.1])

%% local vs. global maximum using max

% create signal
perfdat = -40:.01:120;
signal = 3*sin(perfdat)./perfdat + sin(perfdat-20)./(perfdat-20);

figure(7), clf
plot(perfdat,signal,'k.-','linew',2)
xlabel('X'), ylabel('Y')
title('Two sincs, in synch')

% global maximum
[maxGval,maxGidx] = max(signal);
hold on
plot(perfdat(maxGidx),maxGval,'kv','linew',3,'markersize',15,'markerfacecolor','r')

% local maximum in specified range
range4max = [10 40];

rangeidx = dsearchn(perfdat',range4max');

[maxLval,maxLidx] = max(signal(rangeidx(1):rangeidx(2)));
% is this accurate??
plot(perfdat(maxLidx+rangeidx(1)),maxLval,'kv','linew',3,'markersize',15,'markerfacecolor','g')
zoom on

%% local minima/maxima

% find local maxima
peeks = find(diff(sign(diff(signal)))<0)+1;

[~,idx] = sort(signal(peeks));
peeks(idx(end-1:end)) = [];

plot(perfdat(peeks),signal(peeks),'ro','linew',2,'markersize',10,'markerfacecolor','y')

%% end.

