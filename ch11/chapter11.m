%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 11
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% sine waves

srate = 1000; % sampling rate in Hz
time  = 0:1/srate:5; % units of seconds
f     = 4; % units of Hz
a     = 2; % arbitrary units
th    = pi/2; % in radians

sinewave = a * sin( 2*pi*f*time + th );

figure(1), clf
plot(time,sinewave)
% I think the plot is easier to interpret when the line 
% doesn't touch the axis boundary.
set(gca,'ylim',[min(sinewave) max(sinewave)]*1.05)
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'A sine wave at ' num2str(f) ' Hz' ])

%% complex numbers

x = 4+2i;
x = complex(4,2);

figure(2), clf

% plot the complex number as a vector
plot([0 real(x)],[0 imag(x)],'k-o','linew',3)
set(gca,'xlim',[-5 5],'ylim',[-5 5])

% make the plot look a bit nicer
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
axis square
xlabel('Real axis'), ylabel('Imaginary axis')
title('A complex number as a vector in complex space')

% this text object doesn't really look very nice. Can you improve it?
txt_h = text(real(x),imag(x),sprintf('(%g,%g)',real(x),imag(x)));
set(txt_h,'horizontalAlignment','center')

%% complex sine waves

% complex sine wave
csw = a * exp( 1i*2*pi*f*time + th );

% let's see it in 3D
figure(3), clf
plot3(time,real(csw),imag(csw))

% plot adjustments, etc.
xlabel('Time (s)'), ylabel('real part'), zlabel('imaginary part')
title([ 'A complex sine wave at ' num2str(f) ' Hz' ])
axis square

rotate3d on % now you can click-and-spin the plot

% set viewpoint
view([-55 20])

%% extracting information from the complex dot product

% vector of random numbers
sinewaveX = fft(sinewave)/length(time);
randvect = randn(size(sinewaveX));

% complex dot product 
% The result is complex even though one of the vectors is real-valued!
cdp = sum( randvect.*sinewaveX );

% extract the length of the line (do these lines produce different results?)
linemag = sqrt( real(cdp).^2 + imag(cdp).^2 );
linemag = abs( cdp );
linemag = sqrt( cdp.*conj(cdp) );

% extract angle with respect to the positive real axis (in radians)
% these three lines of code are not always equivalent!
lineangle = atan( imag(cdp) / real(cdp) );
lineangle = atan2(imag(cdp),real(cdp));
lineangle = angle( cdp );

% extract projection onto real axis (same thing as the real part)
realpart = real( cdp );


figure(4), clf

% plot the complex vector
plot([0 real(cdp)],[0 imag(cdp)],'k-o','linew',3)
set(gca,'xlim',[-linemag linemag],'ylim',[-linemag linemag])
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
axis square
xlabel('Real axis'), ylabel('Imaginary axis')
title('A complex number as a vector in complex space')

%% time and frequency domain representations of the same signal

% define signal properties
srate = 200; % sampling rate, sometimes abbreviated Fs
time  = 0:1/srate:5;

% define sine wave properties (hint: try changing these and 
% see what happens to the plot)
frequency = 3;
amplitude = 2;
phase     = pi/3;

% create sine wave
sinewave = amplitude * sin( 2*pi*frequency*time + phase );

figure(5), clf
subplot(211)
plot(time,sinewave)
set(gca,'ylim',[min(sinewave) max(sinewave)]*1.05)
xlabel('Time (s)'), ylabel('Amplitude')
title('TIME DOMAIN')

% the code below uses Fourier transform functions that you'll learn about
% later in this chapter. 

sinewaveX = fft(sinewave)/length(time);
hz = linspace(0,srate/2,floor(length(time)/2)+1);
subplot(212)
bar(hz,2*abs(sinewaveX(1:length(hz))))
set(gca,'xlim',[0 frequency*3],'ylim',[0 max(2*abs(sinewaveX))*1.05])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('FREQUENCY DOMAIN')

%% Many signals are easier to decompose in the frequency domain

% sum several sine waves
sinewave = 10 * sin( 2*pi*5*time  + phase ) + ...
           20 * sin( 2*pi*3*time  + phase ) + ...
            5 * sin( 2*pi*15*time + phase ) + ...
           15 * sin( 2*pi*7*time  + phase );

% plot in the time domain
figure(6), clf
subplot(211)
plot(time,sinewave)
set(gca,'ylim',[min(sinewave) max(sinewave)]*1.05)
xlabel('Time (s)'), ylabel('Amplitude')
title('TIME DOMAIN')


% the code below uses Fourier transform functions that you'll learn about
% later in this chapter. 

sinewaveX = fft(sinewave)/length(time);
hz = linspace(0,srate/2,floor(length(time)/2)+1);
subplot(212)
bar(hz,2*abs(sinewaveX(1:length(hz))))
set(gca,'xlim',[0 20],'ylim',[0 max(2*abs(sinewaveX))*1.05])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('FREQUENCY DOMAIN')

%% The slow discrete Fourier transform

% redefine signal (variables taken from first cell)
signal      = a * sin( 2*pi*f*time + th );

% Fourier transform variables
N           = length(signal); % length of sequence
fourierTime = ((1:N)-1)/N;   % "time" vector
fourierTime = (0:N-1)/N;    % same as previous line!
nyquist     = srate/2;     % Nyquist frequency 

% initialize output matrix of Fourier coefficients
signalX = zeros(size(signal)); 


% These are the actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series (plus DC). 
hz = linspace(0,nyquist,floor(N/2)+1);


% loop over frequencies
for fi=1:N
    
    % create sine wave for this frequency
    fourierSine = exp( -1i*2*pi*(fi-1).*fourierTime );
    
    % compute dot product as sum of point-wise elements
    signalX(fi) = sum( fourierSine.*signal );
    
    % note: this can also be expressed as a vector-vector product
    %fourierCoefs(fi) = fourierSine*signal';
end

% scale Fourier coefficients to original scale
signalX = signalX / N;


% plot time-domain signal
figure(7), clf
subplot(221)
plot(time,signal)
xlabel('Time (s)'), ylabel('Amplitude')
title('Data')

% plot an example of one sine wave from the Fourier transform
subplot(222)
plot(fourierTime,real(exp( -2*pi*1i*(10).*fourierTime )))
xlabel('time (a.u.)'), ylabel('Amplitude')
title('Sine wave')

% plot the amplitude spectrum (magnitude of Fourier coefficients)   
subplot(212)
bar(hz,abs(signalX(1:length(hz)))*2)
set(gca,'xlim',[0 f*3],'ylim',[0 max(abs(signalX))*2.3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% getting close to the Nyquist

% subsample a high-sampling rate sine wave (pretend it's a continuous wave)
srate = 1000;
t = 0:1/srate:1;
f = 10; % frequency of the sine wave Hz
d = sin(2*pi*f*t);


% "Measurement" sampling rates
srates = [15 22 50 200]; % in Hz

figure(8), clf
for si=1:4
    subplot(2,2,si)
    
    % plot 'continuous' sine wave
    plot(t,d), hold on
    
    % plot sampled sine wave
    samples = round(1:1000/srates(si):length(t));
    plot(t(samples),d(samples),'ro-','linew',2)
    
    title([ 'Sampled at ' num2str(srates(si)/f) ' times' ])
    set(gca,'ylim',[-1.1 1.1],'xtick',0:.25:1)
end

%% The fast Fourier transform

% super-easy!
signalX1 = fft( signal )/N;

% plot signalX1 and signalX on top of each other to show that they are identical.

%% more on Fourier coefficients as complex numbers

srate = 1000;
time  = 0:1/srate:1;
freq  = 6;

% create sine waves that differ in power and phase
cos1 = 3 * cos(2*pi*freq*time + 0 );
cos2 = 2 * cos(2*pi*freq*time + pi/6 );
cos3 = 1 * cos(2*pi*freq*time + pi/3 );

% compute Fourier coefficients
fCoefs1 = fft(cos1) / length(time);
fCoefs2 = fft(cos2) / length(time);
fCoefs3 = fft(cos3) / length(time);

hz = linspace(0,srate/2,floor(length(time)/2)+1);

% find the frequency of our sine wave
hz6 = dsearchn(hz',freq);

% let's look at the coefficients for this frequency
disp([ '6 Hz Fourier coefficient for cos1: ' num2str(fCoefs1(hz6)) ])
disp([ '6 Hz Fourier coefficient for cos2: ' num2str(fCoefs2(hz6)) ])
disp([ '6 Hz Fourier coefficient for cos3: ' num2str(fCoefs3(hz6)) ])

%% complex numbers as vectors in a polar plot

% make polar plots of fourier coefficients
figure(9), clf
h(1) = polar([0 angle(fCoefs1(hz6))],[0 2*abs(fCoefs1(hz6))],'r');
hold on
h(2) = polar([0 angle(fCoefs2(hz6))],[0 2*abs(fCoefs2(hz6))],'b');
h(3) = polar([0 angle(fCoefs3(hz6))],[0 2*abs(fCoefs3(hz6))],'k');
set(h,'linewidth',5)
legend({'cos1';'cos2';'cos3'})

%% phase and power information can be extracted via Euler's formula

% extract amplitude using Pythagorian theorem
amp1 = sqrt( imag(fCoefs1).^2 + real(fCoefs1).^2 );
amp2 = sqrt( imag(fCoefs2).^2 + real(fCoefs2).^2 );
amp3 = sqrt( imag(fCoefs3).^2 + real(fCoefs3).^2 );

% extract amplitude using the Matlab function abs
% amp1 = abs( fCoefs1 );
% amp2 = abs( fCoefs2 );
% amp3 = abs( fCoefs3 );

% and again the same thing using the complex conjugate
% amp1 = sqrt( fCoefs1.*conj(fCoefs1) );
% amp2 = sqrt( fCoefs2.*conj(fCoefs2) );
% amp3 = sqrt( fCoefs3.*conj(fCoefs3) );

figure(10), clf

% plot amplitude spectrum
subplot(211)
plot(hz,2*amp1(1:length(hz)),'ro-','linew',3), hold on
plot(hz,2*amp2(1:length(hz)),'bp-','linew',3)
plot(hz,2*amp3(1:length(hz)),'ks-','linew',3)
set(gca,'xlim',[freq-3 freq+3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% and now for phase...

% extract phase angles using trigonometry
phs1 = atan2( imag(fCoefs1) , real(fCoefs1) );
phs2 = atan2( imag(fCoefs2) , real(fCoefs2) );
phs3 = atan2( imag(fCoefs3) , real(fCoefs3) );

% extract phase angles using Matlab function angle
% phs1 = angle( fCoefs1 );
% phs2 = angle( fCoefs2 );
% phs3 = angle( fCoefs3 );


% plot phase spectrum
subplot(212)
plot(hz,phs1(1:length(hz)),'ro-','linew',3), hold on
plot(hz,phs2(1:length(hz)),'bp-','linew',3)
plot(hz,phs3(1:length(hz)),'ks-','linew',3)
set(gca,'xlim',[freq-3 freq+3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% DC reflects the mean offset

% a very simple signal and its fft
signal  = [ 1 0 1 3 -3 1 2 0 -3 0 -1 2 1 -1];
N = length(signal);
signalX = fft(signal) / N;

% get vector of frequencies in Hz
srate = 100; % in Hz, I'm just making up this number
hz = linspace( 0, srate/2, floor(N/2)+1 );

figure(11), clf

% compute the FFT of the same signal with different DC offsets
signalX1 = fft(signal) / N;
signalX2 = fft(signal-mean(signal)) / N;
signalX3 = fft(signal+1) / N;


% plot signals in the time domain
subplot(211)
plot(ifft(signalX1)*N,'bo-'), hold on
plot(ifft(signalX2)*N,'rd-')
plot(ifft(signalX3)*N,'k*-')
xlabel('Tme (ms)'), ylabel('Amplitude')
set(gca,'xlim',[.5 N+.5])


% plot signals in the frequency domain
subplot(212)
plot(hz,2*abs(signalX1(1:length(hz))),'bo-'), hold on
plot(hz,2*abs(signalX2(1:length(hz))),'rd-')
plot(hz,2*abs(signalX3(1:length(hz))),'k*-')
set(gca,'xlim',[-.5 51])

xlabel('Frequencies (Hz)'), ylabel('Amplitude')
legend({'original';'de-meaned';'increased mean'})

%% zero-padding

figure(12), clf

signal  = [ 1 0 1 2 -3 1 2 0 -3 0 -1 2 1 -1];
N = length(signal);

% compute the FFT of the same signal with different DC offsets
signalX1 = fft(signal,N)     / N;
signalX2 = fft(signal,N+10)  / N;
signalX3 = fft(signal,N+100) / N;

% define frequencies vector
frex1 = linspace( 0, .5, floor(length(signalX1)/2)+1 );
frex2 = linspace( 0, .5, floor(length(signalX2)/2)+1 );
frex3 = linspace( 0, .5, floor(length(signalX3)/2)+1 );


% plot signals in the time domain
subplot(211)
plot(ifft(signalX1)*N,'bo-'), hold on
plot(ifft(signalX2)*N,'rd-')
plot(ifft(signalX3)*N,'k*-')
xlabel('Time points (arb. units)')

% plot signals in the frequency domain
subplot(212)
plot(frex1,2*abs(signalX1(1:length(frex1))),'bo-'), hold on
plot(frex2,2*abs(signalX2(1:length(frex2))),'rd-')
plot(frex3,2*abs(signalX3(1:length(frex3))),'k*-')

xlabel('Normalized frequency units'), ylabel('Amplitude')
legend({'"Native" N';'N+10';'N+100'})

%% 2D FFT

pic = imread('saturn.png');

% shrink the image down to 2D
pic = squeeze(mean(pic,3));

% Matlab's 2D FFT function
picX = fft2(pic);

% same result as above
picX1 = fft(pic);
picX1 = fft(picX1,[],2);

figure(13), clf

% show original ("space domain") image
subplot(131)
imagesc(pic)
axis off, axis square
title('''Space'' domain')

% log-power spectrum
subplot(132)
imagesc(log10(abs(picX)))
set(gca,'clim',[2.5 5])
axis off, axis square
title('Original power spectrum')

% same log-power spectrum but shifted
subplot(133)
imagesc(fftshift(log10(abs(picX))))
set(gca,'clim',[2.5 5])
axis off, axis square
title('shifted power spectrum')

%% end.
