%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 20
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% running-mean filter

srate = 1000;
time  = 0:1/srate:6;
ampl  = interp1(0:5,rand(6,1),time,'spline');
f     = 8; % Hz because time is in seconds
noise = .3*randn(size(time));
signal = ampl .* sin(2*pi*f*time) + noise;

k = 7;
filtsig = zeros(size(signal));

for i=k+1:length(time)-k
    filtsig(i) = mean( signal(i-k:i+k) );
end

figure(1), clf
plot(time,signal,time,filtsig)
legend({'original';'filtered'})
title('Illustration of running-mean filter')

%% running-median filter

srate = 1000;
time  = 0:1/srate:6;
ampl  = interp1(0:5,rand(6,1),time,'spline');
f     = 8; % Hz because time is in seconds
noise = 1000*isprime(1:length(time));
signal = ampl .* sin(2*pi*f*time) + noise;

k = 7;
[filtMean,filtMedian] = deal( zeros(size(signal)) );

for i=k+1:length(time)-k
    filtMean(i) = mean( signal(i-k:i+k) );
    
    sortnums = sort( signal(i-k:i+k) );
    filtMedian(i) = sortnums(k);
end

figure(2), clf

subplot(221)
plot(time,signal)
set(gca,'xlim',[1 2],'ylim',[-1.5 1.5])
title('Signal')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(222)
plot(time,filtMean)
title('Mean-filtered')
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'xlim',[1 2],'ylim',[-50 350])


subplot(223)
plot(time,filtMedian)
title('Median-filtered')
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'xlim',[1 2],'ylim',[-1.5 1.5])

zoom on

%% edge in the frequency domain and its time-domain representation

% sharp point (a.k.a. impulse response function, a.k.a. delta function)
N = 400;
X = zeros(N,1);
X(round(N*.05)) = 1;

figure(3), clf
subplot(221)
plot(X)
set(gca,'ylim',[-.01 1.05])
xlabel('Frequency (a.u.)'), ylabel('amplitude')

subplot(223)
plot(real(ifft(X))*length(X))
set(gca,'ylim',[-1.1 1.1])
xlabel('Time (a.u.)'), ylabel('amplitude')


% box-car
X = zeros(N,1);
X(10:30) = 1;

subplot(222)
plot(X)
set(gca,'ylim',[-.01 1.05])
xlabel('Frequency (a.u.)'), ylabel('amplitude')

subplot(224)
plot(real(ifft(X))*length(X))
set(gca,'ylim',[-15 15])
xlabel('Time (a.u.)'), ylabel('amplitude')


%% empirical FWHM of a Gaussian

% create a gaussian
x = -4:.1:4;
gaus = exp(-(x.^2));
gaus = gaus./max(gaus);

% find index of peak
[~,pidx] = max(gaus);
prepeak  = dsearchn(gaus(1:pidx)',.5);
postpeak = pidx-1+dsearchn(gaus(pidx:end)',.5);

% plot gaussian
figure(4), clf
plot(x,gaus,'k-o')
hold on

% plot the empirical points closest to 50% amplitude
plot(x(prepeak),gaus(prepeak),'ro')
plot(x(postpeak),gaus(postpeak),'ro')
plot([ x(prepeak) x(postpeak)],[gaus(prepeak) gaus(postpeak)],'r--')
plot(get(gca,'xlim'),[.5 .5],'k:')
ylabel('Normalized amplitude')

% so what is the FWHM in Hz?

%% create a frequency-domain Gaussian with specified FWHM in hz

fwhm     =  5; % in Hz
centfreq = 14; % center frequency, also in Hz


% set parameters
srate = 1000;
N = 4000;

hz = linspace(0,srate,N);  % frequencies
s  = fwhm*(2*pi-1)/(4*pi); % normalized width 
x  = hz-centfreq;          % shifted frequencies 
gx = exp(-.5*(x/s).^2);    % gaussian 
gx = gx./max(gx);          % gain-normalized 



% compute empirical frequency and standard deviation
idx = dsearchn(hz',centfreq);
empVals(1) = hz(idx);

% find values closest to .5 after MINUS before the peak
empVals(2) = hz(idx-1+dsearchn(gx(idx:end)',.5)) - hz(dsearchn(gx(1:idx)',.5));


% plot
figure(5), clf
plot(hz,gx,'o-')
hold on
plot([hz(dsearchn(gx(1:idx)',.5)) hz(idx-1+dsearchn(gx(idx:end)',.5))],[gx(dsearchn(gx(1:idx)',.5)) gx(idx-1+dsearchn(gx(idx:end)',.5))],'k--')
set(gca,'xlim',[max(centfreq-fwhm*2,0) centfreq+fwhm*2],'ylim',[-.01 1.01]);

title([ 'Requested: ' num2str(centfreq) ', ' num2str(fwhm) ' Hz; Empirical: ' num2str(empVals(1)) ', ' num2str(empVals(2)) ' Hz' ])
xlabel('Frequency (Hz)'), ylabel('Amplitude gain')

%% now we apply the filter to a chirp

% re-create chirp signal
% same length as Gaussian created in previous cell!
sigtime = 0:1/srate:4-1/srate;
n = length(sigtime);

f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);


fwhm     =  3; % in Hz
centfreq =  8;


N = length(sigtime);

hz = linspace(0,srate,N);  % frequencies
s  = fwhm*(2*pi-1)/(4*pi); % normalized width 
x  = hz-centfreq;          % shifted frequencies 
gx = exp(-.5*(x/s).^2);    % gaussian 
gx = gx./max(gx);          % gain-normalized 


% apply the filter
filtsig = 2*real( ifft( fft(signal).*gx ) );

% and plot
figure(6), clf
h=plot(sigtime,signal);
hold on
plot(sigtime,filtsig,'k','linew',2)

%% setting up an FIR filter

% sampling rate
srate = 500;

% shape of the frequency response
fshape = [ 0 1 1 0 ];
frex = [ 0 12 18 srate/2 ]./srate/2;

% plot it.
figure(7), clf

subplot(211)
plot(frex,fshape,'k-o','linew',3,'markersize',8)
xlabel('Frequency (fraction of Nyquist)'), ylabel('Gain')
set(gca,'ylim',[-.02 1.02])

% That doesn't look very good. Actually, it looks pretty awful.
% Let's add some more parameters...

tz = .15; % transition zone, in percent
fbnd = [ 12 18 ]; % frequency boundaries
fshape = [ 0 0 1 1 0 0 ]; % note the difference between this and the previous shape
frex = [ 0 fbnd(1)*(1-tz) fbnd fbnd(2)*(1+tz) srate/2 ];
frex = frex./(srate/2); % norm. to Nyquist

% plot the improved filter response
subplot(212)
plot(frex,fshape,'k-o','linew',3,'markersize',8)
xlabel('Frequency (Hz)'), ylabel('Gain')
set(gca,'ylim',[-.02 1.02])

%% apply filter to linear chirp

% this time we'll use Matlab's chirp function
srate  = 500;
time   = 0:1/srate:15;
frange = [5 20];
x      = chirp(time,frange(1),time(end),frange(2)); % for octave: 'chirp' is in the signal package

% setup the FIR filter again
tz     = .15; % transition zone, in percent
fbnd   = [ 13 15 ]; % freq boundaries
fshape = [ 0 0 1 1 0 0 ];
frex   = [ 0 fbnd(1)*(1-tz) fbnd fbnd(2)*(1+tz) srate/2 ];
frex   = frex./(srate/2); % norm. to Nyquist

% filter order
orderfact = 3;
filtord = round( (orderfact*1000/fbnd(1)) / (srate/1000) );

% create the filter kernel using Matlab's firls function
filtkernel = firls(filtord,frex,fshape); % for octave: also in signal package

% filter the data
fdata = filtfilt(filtkernel,1,x); % for octave: also in signal package


figure(8), clf
subplot(221)
kernelTime = linspace(1,orderfact*1000/fbnd(1),filtord+1);
plot(kernelTime,filtkernel)
set(gca,'xlim',kernelTime([1 end]))
xlabel('Time (ms)')
title('Filter kernel in the time domain')


subplot(222)
hz = linspace(0,srate/2,floor((filtord+1)/2)+1);
powr = abs(fft(filtkernel)).^2;
plot(hz,powr(1:length(hz)))
hold on
plot(frex*srate/2,fshape,'ro-')
set(gca,'xlim',[0 fbnd(end)*3],'ylim',[-.05 1.1])
legend({'Empirical filter';'''Ideal'' filter'})
xlabel('Frequency (Hz)')
title('Filter kernel in the frequency domain')


subplot(212)
plot(time,x,'m'), hold on
plot(time,fdata,'k','linew',2)
% set(gca,'xlim',[6 9])
pan on

%% hilbert transform

% a random signal and its Fourier coefficients
n = 20;
r = randn(n,1);
rx = fft(r);

% a copy that is multiplied by the complex operator
rxi = rx*1i; 

% find indices of positive and negative frequencies
posF = 2:floor(n/2)+mod(n,2); 
negF = ceil(n/2)+1+~mod(n,2):n; 

% rotate Fourier coefficients 
rx(posF) = rx(posF) + -1i*rxi(posF); 
rx(negF) = rx(negF) +  1i*rxi(negF); 

% take inverse FFT 
hilbert_r = ifft(rx);

figure(9), clf

subplot(311)
plot(r,'o-','markersize',9)
hold on
plot(real(hilbert_r),'k*')
plot(real( hilbert(r) ),'gd','markersize',20) % for octave: hilbert function is in the signal package
legend({'original';'real part of Hilbert';'real part of Matlab hilbert'})
title('Real part')


subplot(312)
plot(abs(hilbert_r).^2,'k*-')
hold on
plot(abs(hilbert(r)).^2,'gd','markersize',20)
xlabel('Frequency (a.u.)'), ylabel('Amplitude (a.u.)')
legend({'power from Hilbert';'power from Matlab hilbert'})
title('Power')


subplot(313)
plot(angle(hilbert_r),'k*-')
hold on
plot(angle(hilbert(r)),'gd','markersize',20)
xlabel('Frequency (a.u.)'), ylabel('Phase (rad.)')
legend({'phases from Hilbert';'phases from Matlab hilbert'})
title('Phases')

%% end.
