%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 23
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% import data

% import data from previous chapter
fid = fopen('../ch22/03-04-11-aa_sig1_spikes.dat');

spikenum = 1;
trialnum = 0;
spikesparse = zeros(100,2);
spikesfull  = zeros(1,1001);

while ~feof(fid)
    
    % get a line of code
    d=fgetl(fid);
    
    % if it's an empty line, skip to the next line
    if isempty(d), continue; end
    
    % if this is an inter-trial comment, increment trial number 
    % and then skip forward
    if d(1)=='#'
        if strcmpi(d(1:5),'# int')
            trialnum = trialnum+1;
        end
        continue;
        
        
    % otherwise, it's real data; grab the spike time 
    % and update the matrices
    else
        spiketime = ceil(sscanf(d,'%g')/1000);
        
        spikesparse(spikenum,:) = [trialnum spiketime];
        spikesfull(trialnum,spiketime) = 1;
        spikenum = spikenum+1;
    end
end
fclose(fid);


figure(1), clf
plot(spikesparse(:,2),spikesparse(:,1),'m.')
xlabel('Time (ms)'), ylabel('Trials')

%% spike timing with full matrices

win = 50; % in ms and also in indices (only in this dataset!)
spikerhyth = zeros(1,win*2+1);
n = 0;

for triali=1:trialnum
    
    spikes = find(spikesfull(triali,:));
    % need an extra line of code here...
    
    for si=1:length(spikes)
        spikerhyth = spikerhyth + spikesfull(triali,spikes(si)-win:spikes(si)+win);
        n = n+1;
    end % end spike loop
end % end trial loop

% divide by N to finalize average
spikerhyth = spikerhyth./n;

figure(2), clf
plot(-win:win,spikerhyth,'rs-')
set(gca,'xlim',[-10 10])

%% spike timing with sparse matrices

win = 50; % in ms and also in indices
spikerhyth = zeros(1,win*2+1);
n = 0;

for triali=1:max(spikesparse(:,1))
    
    spikes = spikesparse(spikesparse(:,1)==triali,2);
    
    for si=1:length(spikes)
        tempspikes = spikes-spikes(si);
        tempspikes(tempspikes<-win | tempspikes>win) = [];
        
        spikerhyth(tempspikes+win+1) = spikerhyth(tempspikes+win+1)+1;
        
        n = n+1;
    end % end spike loop
end % end trial loop
spikerhyth = spikerhyth./n;


% plot on top of the previous result to see how they compare
figure(2), hold on
plot(-win:win,spikerhyth,'ko-')
set(gca,'xlim',[-10 10])
legend({'From full matrix';'From sparse matrix'})

%% using autocorrelation in the frequency domain

spikerhyth = zeros(1,win*2+1);

for triali=1:trialnum
    
    % mean subtract data from this trial
    meansub = bsxfun(@minus,spikesfull(triali,:),mean(spikesfull(triali,:)));
    
    % take FFT and iFFT of power spectrum
    temp = fft(meansub);
    temp = ifft( abs(temp).^2 );
    
    % chop out parts we want
    temp = temp([end-win+1:end 1:win+1]);
    
    % sum autocorrelation
    spikerhyth = spikerhyth + temp./max(temp);
end

% again, divide by N
spikerhyth = spikerhyth./triali;


% and plot
figure(3), clf
plot(-win:win,spikerhyth,'rd-','markersize',10)
set(gca,'xlim',[-10 10])

%% the previous cell can also be done more efficiently

% Using one really long line of code! So I guess it depends on your definition of "efficient."
spikerhyth = mean( ifft( abs(fft( bsxfun(@minus,spikesfull,mean(spikesfull,2)),[],2 )).^2,[],2 ) ,1) / trialnum;

figure(3), hold on
plot(-win:win,spikerhyth([end-win+1:end 1:win+1]),'m')
set(gca,'xlim',[-10 10])

%% cross-neuron spike timing correlation

% try running datasets 115 and 111
load times_090425blk10_ch115.mat

% spike timing with sparse matrices
win = 50; % in ms
spikerhyth = zeros(1,win*2+1);
n = 0;

spikesClus1 = cluster_class(cluster_class(:,1)==1,2);
spikesClus2 = cluster_class(cluster_class(:,1)==2,2);

for si=1:length(spikesClus1)
    tempspikes = round(spikesClus2-spikesClus1(si));
    tempspikes(tempspikes<-win | tempspikes>win) = [];
    
    spikerhyth(tempspikes+win+1) = spikerhyth(tempspikes+win+1)+1;
    
    n = n+1;
end % end spike loop
spikerhyth = spikerhyth./n;

figure(4), clf
plot(-win:win,spikerhyth,'kd-','markersize',10)
set(gca,'xlim',[-10 10])
xlabel('Relative spike times (ms)')
ylabel('Probability')

%%



%% hc-5 for spike-lfp coherence

% data are taken from:
% https://crcns.org/data-sets/hc/hc-5/about-hc-5/
load spikefieldData

wherespikes = find(spikeTimes);
win = 200; % indices, not ms!
wherespikes(wherespikes<win | wherespikes>length(lfp)-win) = [];
spikeLFP = zeros(length(wherespikes),win*2+1);

for si=1:length(wherespikes)
    spikeLFP(si,:) = lfp(wherespikes(si)-win:wherespikes(si)+win);
end


figure(5), clf
timevec = (-win:win)/(1000/srate);
plot(timevec,mean(spikeLFP))
hold on
plot([0 0],get(gca,'ylim'),'k')
set(gca,'xlim',timevec([1 end]))
xlabel('Peri-spike time (ms)'), ylabel('Voltage (mV)')

%
figure(6), clf
subplot(121)
h = plot(timevec,spikeLFP(1:100:end,:));
hold on
plot(timevec,mean(spikeLFP),'k','linew',4)
set(gca,'xlim',timevec([1 end]))
set(h,'color',[.7 .7 .7])
xlabel('Peri-spike time (ms)'), ylabel('Trials')


subplot(122)
imagesc(timevec,[],spikeLFP)
set(gca,'clim',[-500 500],'xlim',timevec([1 end]))
xlabel('Peri-spike time (ms)'), ylabel('Trials')

%% same but filtered

frex = 8; % hz

wavetime = -2:1/srate:2;
halfwave = floor(length(wavetime)/2)+1;


nData = length(lfp);
nWave = length(wavetime);
nConv = nData + nWave - 1;

% create wavelet
nCycl = 6;
gausS = nCycl / (2*pi*frex);

wavelet  = exp(2*1i*pi*frex*wavetime + (-wavetime.^2)/(2*gausS^2) );
waveletX = fft(wavelet,nConv);
waveletX = waveletX./max(waveletX);

%% convolution

dataX = fft(lfp,nConv);
as = ifft( dataX .* waveletX );

% cut off edges
as = as(halfwave-1:end-halfwave);


% now
sfc = abs(mean(exp(1i*angle(as(logical(spikeTimes))))));


wherespikes = find(spikeTimes);
win = 200; % indices, not ms!
wherespikes(wherespikes<win | wherespikes>length(as)-win) = [];
spikeLFP = zeros(length(wherespikes),win*2+1);

for si=1:length(wherespikes)
    spikeLFP(si,:) = real(as(wherespikes(si)-win:wherespikes(si)+win));
end

figure(7), clf
timevec = (-win:win)/(1000/srate);
plot(timevec,mean(spikeLFP))
set(gca,'xlim',timevec([1 end]))

%
figure(8), clf
subplot(121)
h = plot(timevec,spikeLFP(1:100:end,:));
hold on
plot(timevec,mean(spikeLFP),'k','linew',4)
set(gca,'xlim',timevec([1 end]))
set(h,'color',[.7 .7 .7])


subplot(122)
imagesc(timevec,[],spikeLFP)
set(gca,'clim',[-200 200])

%% now smooth the image

gx = -10:10;

gaus2d = zeros(length(gx));

s = 10;

for xi=1:length(gx)
    for yi=1:length(gx)
        gaus2d(xi,yi) = exp( -( gx(xi)^2+gx(yi)^2 ) / (2*s^2) );
    end
end


figure(9), clf

subplot(211)
imagesc(gaus2d)
axis image

subplot(212)
smo = conv2(spikeLFP,gaus2d,'same');
imagesc(timevec,[],smo./sum(gaus2d(:)))
set(gca,'clim',[-100 100])

%%


