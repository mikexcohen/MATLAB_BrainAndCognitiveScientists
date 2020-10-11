%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 19
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% complex Morlet wavelet

srate = 1000;
wavtime = -2:1/srate:2;
frex = 6.5;
s = 5/(2*pi*frex);
csine = exp(2*1i*pi*frex*wavtime);
gaus = exp( -(wavtime.^2) / (2*s^2) );
cmw = csine .* gaus;

figure(1), clf
plot3(wavtime,real(cmw),imag(cmw))
xlabel('Time'), ylabel('Real part'), zlabel('Imaginary part')
axis image
title('Complex Morlet wavelet in the time domain')
rotate3d on


figure(2), clf
hz = linspace(0,srate/2,floor(length(wavtime)/2)+1);
cmwX = fft(cmw)/length(wavtime);
plot(hz,abs(cmwX(1:length(hz)))*2)
set(gca,'xlim',[0 frex*3])
xlabel('Frequencies (Hz)'), ylabel('Amplitude')
title('Power spectrum of Morlet wavelet')

%% Convolution!

% some signal parameters
sigtime = 0:1/srate:10;
n = length(sigtime);

% create a linear chirp
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);

% convolution parameters
nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

% the meat and potatoes of convolution
as = ifft( fft(cmw,nConv) .* fft(signal,nConv) );
as = as(nHfkn:end-nHfkn+1);


% plot the time-domain signal (chirp)
figure(3), clf
subplot(211)
plot(sigtime,signal)
set(gca,'ylim',[-1.1 1.1])
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal')

% and the power time series at one frequency
subplot(212)
plot(sigtime,abs(as).^2)
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Time-frequency power (' num2str(frex) ' Hz)' ])

%% now a time-frequency plane

% new frequency parameters
nFrex = 30;
frex = linspace(1,15,nFrex);

% width of Gaussian
s = linspace(4,12,nFrex) ./ (2*pi.*frex);

% initialize output matrix
tf = zeros(nFrex,length(signal));

for fi=1:nFrex
    
    % Take the FFT of the signal.
    % Does it make sense to have this code *inside* the frequency loop?
    sigX = fft(signal,nConv);
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime + ...
              -(wavtime.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    tf(fi,:) = abs( as(nHfkn:end-nHfkn+1) )*2;
end

% and plot
figure(4), clf
contourf(sigtime,frex,tf,40,'linecolor','none')

%% check the frequency representation of the wavelets

figure(5), clf

% hz vector
hz = linspace(0,srate/2,floor(length(wavtime)/2)-1);

for fi=1:nFrex
    
    % create wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    cmwX = fft(cmw,nConv);
    
    % plot wavelet and its power spectrum
    subplot(211), plot(real(cmw))
    subplot(212), plot(hz,2*abs(cmwX(1:length(hz))));
    
    title([ 'Frequency = ' num2str(frex(fi)) ])
    set(gca,'xlim',[0 80])
    pause % only the user can continue
end

%% try with real EEG data, and illustrate the timetrial

% load EEG data
load ../sampleEEGdata.mat

% create supertrial. What is the size of variable data
% in the next two lines? How could you recover the time-by-trials matrix?
data = squeeze(EEG.data(47,:,:));
data = reshape(data,1,[]);

% create wavelet
wavtime = -2:1/EEG.srate:2;
frex    = 6.5;
s       = 5/(2*pi*frex); % what is the number of cycles here?
csine   = exp(2*1i*pi*frex*wavtime);
gaus    = exp( -(wavtime.^2) / (2*s^2) );
cmw     = csine .* gaus;


% step 1: define convolution parameters
nData = length(data); % or: EEG.pnts*EEG.trials
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

% step 2: take FFTs
dataX = fft(data,nConv);
cmwX  = fft(cmw,nConv);

% step 3: normalize kernel
cmwX  = cmwX./max(cmwX); % normalize to get the right units

% step 4: point-wise multiply and take iFFT
as = ifft( dataX.*cmwX );

% step 5: trim wings
as = as(nHfkn:end-nHfkn+1);

% step 5.5: check the size of variable as. Then try plotting it before reshape:
% plot3(1:length(as),real(as),imag(as))

% new step 6: reshape back to time-by-trials
as = reshape(as,EEG.pnts,EEG.trials);


% now plot!
figure(6), clf

subplot(211)
plot(EEG.times,abs(as)*2)
xlabel('Time (ms)'), ylabel('Power')
title([ 'Power at ' num2str(frex) ' Hz from all trials' ])

subplot(212)
plot(EEG.times,mean( abs(as)*2 ,2))
xlabel('Time (ms)'), ylabel('Power')
title([ 'Average power at ' num2str(frex) ' Hz over all trials' ])


% and plot ITPC
figure(7), clf

subplot(211)
plot(EEG.times,angle(as))
xlabel('Time (ms)'), ylabel('Phase angles (rad.)')
title([ 'Phase time series at ' num2str(frex) ' Hz from all trials' ])

subplot(212)
plot(EEG.times,abs( mean(exp(1i*angle(as)),2) ))
xlabel('Time (ms)'), ylabel('Phase clustering')
title([ 'Phase clustering time series at ' num2str(frex) ' Hz over trials' ])

%% edge effects in time-frequency analyses

srate = 1000;
sigtime = 0:1/srate:2;
signal = zeros(size(sigtime));
signal(dsearchn(sigtime',.6):dsearchn(sigtime',1.4)) = 1;

nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

as = ifft( fft(cmw,nConv) .* fft(signal,nConv) );
as = as(nHfkn:end-nHfkn+1);


nFrex = 50;
frex  = linspace(1,srate/10,nFrex);
s     = linspace(4,12,nFrex) ./ (2*pi.*frex);
sigX  = fft(signal,nConv);
tf    = zeros(nFrex,length(signal));

for fi=1:nFrex
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    tf(fi,:) = abs( as(nHfkn:end-nHfkn+1) )*2;
end

figure(8), clf

subplot(211)
plot(sigtime,signal)
set(gca,'ylim',[-.05 1.05])

subplot(212)
contourf(sigtime,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .1])

%% reflection

t = -1:1/srate:3;
g = [ diff( exp(-t.^2/.1) ) 0];
reflex = [ g(end:-1:1) g g(end:-1:1) ];

figure(9),clf

% plot original signal
subplot(211)
plot(g)
set(gca,'xlim',[0 numel(reflex)]-numel(g))

% plot reflected signal
subplot(212)
plot(length(g)+1:2*length(g),g,'ko')
hold on
plot(reflex)
set(gca,'xlim',[0 numel(reflex)])
legend({'original';'reflected'})

%% short-time FFT

% re-create chirp signal
srate = 1000;
sigtime = 0:1/srate:10;
n = length(sigtime);

f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);



% STFFT parameters
fftWidth_ms = 1000; % FFT width in ms

% convert width from ms to indices
fftWidth = round(fftWidth_ms/(1000/srate)/2);
Ntimesteps = 50; % number of time widths
ct = round(linspace(fftWidth+1,n-fftWidth,Ntimesteps));

% figure out Hz vector for our FFT width
hz = linspace(0,srate/2,fftWidth-1);

hwin = .5*(1-cos(2*pi*(1:fftWidth*2)/ (fftWidth*2-1)));

tf = zeros(length(hz),length(ct));
for ti=1:length(ct)
    tdat = signal(ct(ti)-fftWidth:ct(ti)+fftWidth-1);
    x = fft(hwin.*tdat)/fftWidth;
    tf(:,ti) = 2*abs(x(1:length(hz)));
end


figure(9), clf
subplot(211), plot(sigtime,signal)
set(gca,'ylim',[-1.1 1.1])

subplot(212)
contourf(ct,hz,tf,10,'linecolor','none')
set(gca,'ylim',[0 max(f)*1.5],'clim',[0 1])

%% baseline normalization using dB

time = -2:1/100:5; % 100 Hz sampling rate
basetime = [-1.5 -.5];
baseidx = dsearchn(time',basetime');

mothersig = sin(2*pi*5* (time-3) )./(time-3);
mothersig(~isfinite(mothersig)) = max(mothersig);

figure(10), clf
subplot(311)
plot(time,mothersig)
xlabel('Time (s)'), ylabel('Amplitude')
title('"Mother" signal')


% create two children signals by adding different DC offsets
sig1 = mothersig + 2000;
sig2 = mothersig +  200;

% compute the average per signal
basePow1 = mean(sig1(:,baseidx(1):baseidx(2)),2);
basePow2 = mean(sig2(:,baseidx(1):baseidx(2)),2);

% dB-normalization
sig1DB = 10*log10( bsxfun(@rdivide,sig1,basePow1) );
sig2DB = 10*log10( bsxfun(@rdivide,sig2,basePow2) );

% plot the two signals in the same plot. Easy to compare them?
subplot(312)
plot(time,sig1,time,sig2)
xlabel('Time (sec.)'), ylabel('Raw amplitude')

% and show the normalized signals
subplot(313)
plot(time,sig1DB,time,sig2DB)
legend({'sig1';'sig2'})
xlabel('Time (sec.)'), ylabel('dB amplitude')
title('Normalized signals')

%% plotting time-domain response on top of TF plot

load sampleEEGdata.mat

% extract time-frequency power here...

%% produce a contour plot of TF data

figure(11), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'ydir','normal','xlim',[-300 1300],'clim',[-3 3])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

% now you can rescale the ERP (time-domain trial average)
erp = squeeze(mean(EEG.data(23,:,:),3));
erp = (erp-min(erp))./max(erp-min(erp));
yscale = get(gca,'ylim');
erp = erp*(yscale(2)-yscale(1))+yscale(1);

hold on
plot(EEG.times,erp,'k','linew',2)

%% end.
