%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 18
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% frequency resolution

srate = 739;  % sampling rate in Hz
n     = 1234; % number of time points

% compute vector of frequencies
hz = linspace(0,srate/2,floor(n/2)+1);

% frequency resolution, measured in several different ways:
% average derivative of frequencies in Hz
freqres = mean(diff(hz));

% difference between any two frequencies
freqres = hz(2)-hz(1);

% ratio of sampling rate to number of data points
freqres = srate/n;

%% finding the right N given a desired frequency resolution

srate   = 1000;
freqres = .25; % in Hz

% notice that this line is the same as (freqres=srate/n) but solved for n
nFFT    = ceil( srate/freqres );

%% difference between sampling rate and number of time points for Fourier frequencies

% This cell illustrates the difference between 
% frequency resolution and the Nyquist frequency.

% temporal parameters
srates  = [100 100 1000];
timedur = [  1  10    1];

% a few more specifications
freq    =     5; % in Hz
colors  = 'kmb';
symbols = 'op.';


figure(1), clf
legendText = cell(size(timedur));
for parami=1:length(colors)
    
    % define sampling rate in this round
    srate = srates(parami); % in Hz
    
    % define time
    time = -1:1/srate:timedur(parami);
    
    % create signal (Morlet wavelet)
    signal = cos(2*pi*freq.*time) .* exp( (-time.^2) / .05 );
    
    % compute FFT and normalize
    signalX = fft(signal)/length(signal);
    signalX = signalX./max(signalX);
    
    % define vector of frequencies in Hz
    hz = linspace(0,srate/2,floor(length(signal)/2)+1);
    
    
    % plot time-domain signal
    subplot(211)
    plot(time,signal,[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami)), hold on
    set(gca,'xlim',[-1 1])
    xlabel('Time (s)'), ylabel('Amplitude')
    title('Time domain')
    
    % plot frequency-domain signal
    subplot(212), hold on
    plot(hz,abs(signalX(1:length(hz))),[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami))
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
    title('Frequency domain')
    
    legendText{parami} = [ 'srate=' num2str(srates(parami)) ', N=' num2str(timedur(parami)+1) 's' ];
end

legend(legendText)
zoom on

%% edges in the time domain produce "artifacts" in the frequency domain

% create time series with sharp edge
ts = zeros(100,1);
ts(48:52) = 10;


figure(2), clf

% plot time-domain signal
subplot(211)
plot(ts,'k','linew',3)
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
set(gca,'ylim',[-.5 10.5])

% plot power spectrum of signal
subplot(212)
plot(abs(fft(ts)/100).^2,'ro','linew',1.5,'markerface','g')
xlabel('Frequencies (a.u.)'), ylabel('Power (a.u.)')

%% illustration of tapering

figure(3), clf

n = 100;
r = randn(1,n);
hannwin = .5*( 1-cos(2*pi*(0:n-1)/(n-1)) );
plot(r), hold on
plot(r.*hannwin,'r')
plot(hannwin,'k--')

xlabel('Time (a.u.)')
ylabel('Amplitude (a.u.)')
title('Tapering tapers the data')

legend({'original data';'tapered data';'Hann taper'})

%% epoching continuous data

srate = 512; % in hz
n     = 21*srate; % 21 seconds of â€œactivityâ€?

% continuous data
data  = randn(1,n);

% desired length of epoch
epochLms  = 2345; % in ms
% and convert to indices (time units)
epochLidx = round( epochLms / (1000/srate) );

% how many epochs will we have?
nEpochs = floor(n/epochLidx);

% finally, create the epochs using however much data is possible
epochs  = reshape( data(1:nEpochs*epochLidx) ,nEpochs,epochLidx);

%% using the reshape function

% original data
data = 1:12;

% various reshaping options
reshape(data,1,[])
reshape(data,[],1)
reshape(data,3,4)
reshape(data,4,3)
reshape(data,4,4)

%% reshaping problems

nchans = 10;
data = randn(nchans,n);
epochs = reshape(data(:,1:nEpochs*epochLidx),[nchans nEpochs epochLidx]);


figure(4), clf

subplot(211)
plot(data(1,1:100),'r'), hold on
plot(squeeze(epochs(1,1,1:100)))
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
title('Uh oh.')

subplot(212)
epochs1 = reshape(data(:,1:nEpochs*epochLidx),[nchans epochLidx nEpochs]);
plot(data(1,1:100),'r','linew',2), hold on
plot(squeeze(epochs1(1,1:100,1)),'ko--')
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
title('Looks better.')

%% finding the right frequency

srate = 1000;
nFFT  = 5128;
hz    = linspace(0,srate/2,floor(nFFT/2)+1);

% does this give frequencies between 8 Hz and 12 Hz?
hz([8 12])


figure(5), clf

desfreq = 8; % desired frequency in Hz

subplot(131)
plot(hz), hold on
plot(get(gca,'xlim'),[desfreq desfreq],'k--')
xlabel('Indices into variable hz')
ylabel('Actual frequencies in Hertz')
set(gca,'xlim',[0 400],'ylim',[0 30])
title('hz')
axis square


subplot(132)
plot(hz-desfreq), hold on
plot(get(gca,'xlim'),[0 0],'k--')
xlabel('Indices into variable hz')
ylabel([ 'Hertz relative to ' num2str(desfreq) ])
set(gca,'xlim',[0 400],'ylim',[-30 30])
title('hz-desfreq')
axis square

subplot(133)
plot(abs(hz-desfreq)), hold on
plot(get(gca,'xlim'),[0 0],'k--')
xlabel('Indices into variable hz')
ylabel([ 'Hertz relative to ' num2str(desfreq) ])
set(gca,'xlim',[0 400],'ylim',[-30 30])
title('abs(hz-desfreq)')
axis square

%% power spectrum of real EEG data before and after tapering

% load EEG data (one electrode from a human during resting-state)
load EEGrestingState
npnts = size(eegdata,1);

dataX = fft(eegdata,[],1)/npnts;
hz = linspace(0,srate/2,floor(npnts/2)+1);

figure(5), clf
subplot(211)
pow = mean( abs(dataX(1:length(hz),:)).^2 ,2);
plot(hz,pow,'k.-')

subplot(212)
plot(hz,pow./max(pow),'k.-')



% uh oh, no tapering. Better redo:
hannwin = .5 - cos(2*pi*linspace(0,1,npnts))/2;
dataTaper = bsxfun(@times,eegdata,hannwin');
dataX = fft(dataTaper,[],1)/npnts;

subplot(211), hold on
pow = mean( abs(dataX(1:length(hz),:)).^2 ,2);
plot(hz,pow,'r.-')

xlabel('Frequency (Hz)'), ylabel('Power (\muV^2)')
set(gca,'xlim',[0 55])
legend({'Raw time series';'Tapered time series'})

subplot(212), hold on
plot(hz,pow./max(pow),'r.-')

xlabel('Frequency (Hz)'), ylabel('Power (norm.)')
set(gca,'xlim',[0 55])
legend({'Raw time series';'Tapered time series'})

%% extracting power from a specified frequency band

desfreq = [8 12]; % frequency boundaries for averaging
idx = dsearchn(hz',desfreq');
pow = mean( abs(dataX(idx(1):idx(2))).^2 );

%% temporal non-stationarities and the FFT

% define some signal parameters
srate = 1000;
t = 0:1/srate:10;
n = length(t);
f = 3; % frequency in Hz

% create time-increasing amplitude, and stationary amplitude
ampl1 = linspace(1,10,n);
% the next line creates a randomly varying amplitude, FYI
% ampl1 = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
ampl2 = mean(ampl1);

% and create the sine waves
signal1 = ampl1 .* sin(2*pi*f*t);
signal2 = ampl2 .* sin(2*pi*f*t);

% FFTs etc.
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2)+1);


figure(6), clf
subplot(231)
h(1) = plot(t,signal2);
hold on
plot(t,signal1,'k')
xlabel('Time'), ylabel('amplitude')

subplot(234)
h(2) = plot(hz,2*abs(signal2X(1:length(hz))),'o-');
hold on
plot(hz,2*abs(signal1X(1:length(hz))),'ks-','markerface','k')

xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[1 7])
legend({'Stationary';'Non-stationary'})

%% frequency non-stationarity

% create chirp and frequency-stationary signal
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal1 = sin(2*pi.*mean(ff).*t);
signal2 = sin(2*pi.*ff.*t);

% FFTs, etc.
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2));

% and plot
subplot(232)
h(3) = plot(t,signal2);
hold on
plot(t,signal1,'k')
xlabel('Time'), ylabel('amplitude')
set(gca,'ylim',[-1.1 1.1])

subplot(235)
plot(hz,2*abs(signal2X(1:length(hz))),'k.-'), hold on
h(4) = plot(hz,2*abs(signal1X(1:length(hz))),'.-');
xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[0 20])


% Gray lines for a cloudy day. Feel free to adjust 
% the color to match your local weather.
set(h,'color',.6*ones(3,1));

%% even weirder stationarities

% triangle-ish frequency modulation
freqTS = abs(mod(t,2)-1)*10;
meanF = mean(freqTS);
k = 2*pi/srate;
y = sin(2*pi.*meanF.*t + k*cumsum(freqTS-meanF));


subplot(233), plot(t,y)
p = abs(fft(y)/n)*2;
set(gca,'ylim',[-1.1 1.1])

subplot(236), plot(hz,p(1:length(hz)))
set(gca,'xlim',[0 20])

%%

srate = 1000; 
time  = 0:1/srate:5;

% create a signal at 10 Hz with amplitude modulated at 1 Hz
sig = sin(2*pi*time) .* sin(2*pi*10*time);

clf

% plot the time-domain signal
subplot(211)
plot(time,sig)
title('Time domain')
xlabel('Time (s)'), ylabel('amplitude')

% and plot it in the frequency domain
subplot(212)
plot(linspace(0,srate,length(time)), abs(fft(sig)/length(sig)).^2)
set(gca,'xlim',[5 15])
title('Frequency domain')
xlabel('Frequency (Hz)'), ylabel('Power')

% questions: Why don't you see the 10 Hz component in the power spectrum?
%            What happens if you add 1 to the amplitude-modulating function (use parentheses)?

%% spectral coherence

% define parameters
srate = 1000;
t     = 0:1/srate:9;
n     = length(t);
hz    = linspace(0,srate/2,floor(n/2));

% define frequencies
f = [10 14 8];

% define 'base' signals
f_ts1 = (2*pi*cumsum(5*randn(1,n)))/srate;
f_ts2 = (2*pi*cumsum(5*randn(1,n)))/srate;
f_ts3 = (2*pi*cumsum(5*randn(1,n)))/srate;

% notice how the different signals contain different
% parts of the 'base' signals
sigA = sin(2*pi.*f(1).*t + f_ts1) + randn(size(t));
sigB = sigA + sin(2*pi.*f(2).*t + f_ts2);
sigA = sigA + sin(2*pi.*f(3).*t + f_ts3);
% sigB = sigA;

% spectral coherence
sigAx = fft(sigA)/n;
sigBx = fft(sigB)/n;
specX = abs(sigAx.*conj(sigBx)).^2;
% specX = specX./( abs(sigAx).^2 .* abs(sigBx).^2 );

%%% Note: The normalization is not valid when there is only a single trial. 
%%%  The data should be epoched, or you can use the function mscohere, which
%%%  will epoch continuous data (also try mscohere forcing a single window to
%%%  confirm that the normalized coherence values are all 1's).


figure(7), clf

% plot time-domain signals
subplot(231)
plot(t,sigA)
set(gca,'xlim',t([1 end]))
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal A')

subplot(232)
plot(t,sigB)
set(gca,'xlim',t([1 end]))
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal B')


% then plot power spectra
subplot(234)
plot(hz,abs(sigAx(1:length(hz))).*2)
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Signal A')

subplot(235)
plot(hz,abs(sigBx(1:length(hz))).*2)
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Signal B')

% finally, plot coherence
subplot(236)
plot(hz,abs(specX(1:length(hz))))
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Coherence')
title('Spectral coherence: A-B')

%% end.

