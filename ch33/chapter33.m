%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 33
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% create the to-be-animated signal

% signal properties
srate = 1000;
timevec = 0:1/srate:5;

% define the frequency time series and amplitude modulation
freqTS = abs(interp1(linspace(timevec(1),timevec(end),10), 5*rand(1,10),timevec,'spline'));
ampmod = abs(interp1(linspace(timevec(1),timevec(end),10),10*rand(1,10),timevec,'spline'));

% finish defining the frequency-nonstationary signal
centfreq = mean(freqTS);
k = (centfreq/srate)*2*pi/centfreq;
y = ampmod .* sin(2*pi.*centfreq.*timevec + k*cumsum(freqTS-centfreq));

figure(1), clf

% You can view the signal in the "boring" way
plot(timevec,y)

% but I prefer this way:
as = hilbert(y); % octave: hilbert is in the signal package (or write your own function from the book!)
plot3(timevec,real(as),imag(as))
rotate3d on; axis square

%% now let's make a movie

% some movie parameters
vidspeed = 1;
timelag  = 500;

% setup the time course plots
figure(2), clf
subplot(311)
as_h = plot(as,'k','linew',2);
xlabel('Real'), ylabel('Imaginary')
set(gca,'xlim',[min(real(as)) max(real(as))],'ylim',[min(imag(as)) max(imag(as))])
axis square

subplot(312)
bp_h = plot(timevec,real(as),'k','linew',2);
set(gca,'xlim',timevec([1 end]),'ylim',[min(real(as)) max(real(as)) ])
xlabel('Time (ms)'), ylabel('Amplitude')

subplot(313)
pw_h = plot(timevec,abs(as).^2,'k','linew',2);
set(gca,'xlim',timevec([1 end]),'ylim',[min(abs(as).^2) max(abs(as).^2) ])
xlabel('Time (ms)'), ylabel('Power')


for ti=1:vidspeed:length(timevec)
    
    % draw complex values in polar space
    set(as_h,'XData',real(as(max(1,ti-timelag):ti)),'YData',imag(as(max(1,ti-timelag):ti)))
    
    % update cartesian plots
    set(bp_h,'XData',timevec(max(1,ti-timelag):ti),'YData',real(as(max(1,ti-timelag):ti)))
    set(pw_h,'XData',timevec(max(1,ti-timelag):ti),'YData',abs(as(max(1,ti-timelag):ti)).^2)
end

%% gabor movie

% setup size limits of the Gabor
lims = [-31 31];
[x,y] = ndgrid(lims(1):lims(2),lims(1):lims(2));

% number of frames in the movie
nFrames = 20;

% time-varying parameters
phases  = linspace(0,pi,nFrames);
rotates = linspace(0,pi/2,nFrames);
widths  = [ linspace(5,12,nFrames/2) linspace(12,5,nFrames/2) ];


% setup the figure
figure(3), clf
subplot(221)
gaus_h = imagesc(randn(size(x)));
title('Gaussian'), axis off, axis square

subplot(222)
sine_h = imagesc(randn(size(x)));
title('Sine'), axis off, axis square

subplot(212)
gabr_h = imagesc(randn(size(x)));
title('Gabor'), axis off, axis square
set(gca,'clim',[-1 1])


% and... action!
for framei=1:nFrames
    
    % define gaussian width
    width = widths(framei);
    
    % define rotate phase and adjust x,y accordingly
    rotphase = rotates(framei);
    xp = x*cos(rotphase) + y*sin(rotphase);
    yp = y*cos(rotphase) - x*sin(rotphase);
    
    
    % define and display the Gaussian
    gaus2d = exp( -(xp.^2 + yp.^2) / (2*width^2) );
    set(gaus_h,'CData',gaus2d);
    
    % define and display the sine wave
    sine2d = sin(2*pi*.05*xp + phases(framei));
    set(sine_h,'CData',sine2d);
    
    % define and display the gabor patch
    set(gabr_h,'CData',sine2d .* gaus2d);
    
    pause(.1)
end

%% wavelet convolution

% load data
load data4heads

% frequencies in Hz to show
frequencies = [6 11];

% initialize wavelet family
time  = -1:1/EEG.srate:1;
s     = 4.5./(2*pi.*frequencies);
nKern = length(time);
nConv = EEG.pnts*EEG.trials + nKern - 1;
cmwX  = zeros(length(frequencies),nConv);

% Loop through frequencies and make a family of wavelets.
for fi=1:length(frequencies)
    % create complex Morlet wavelets and normalize
    cmw = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*s(fi)^2));
    tmp = fft(cmw,nConv)/length(time);
    cmwX(fi,:) = tmp./max(tmp);
end

% initialize output matrix
tf = zeros(EEG.nbchan,2,EEG.pnts);
basetimeidx = dsearchn(EEG.times',[-500 -100]');

% run convolution
for chani=1:EEG.nbchan
    
    % FFT of EEG data for wavelet convolution
    eegX = fft(reshape(EEG.data(chani,:,:),1,[]),nConv);
    
    for fi=1:length(frequencies)
        
        % wavelet convolution
        ift   = ifft(eegX.*cmwX(fi,:));
        temp  = ift(floor(nKern/2)+1:end-round(nKern/2)+1);
        
        % reshape to 2D, then get trial average (all in one line)
        temp = mean(abs(reshape(temp,EEG.pnts,EEG.trials)).^2,2);
        
        % dB-normalize and put in matrix
        tf(chani,fi,:) = 10*log10(bsxfun(@rdivide,temp,mean(temp(basetimeidx))));
    end
end

%% that weird dream that David Lynch might have had.

figure(4), clf
set(gcf,'color','k','InvertHardcopy','off')

% number of movie frames
nframes = 100;

% create a nice circular colormap
cmap=(1+[cos(linspace(0,pi*2,100)); sin(linspace(0,pi*2,100)); cos(linspace(0,pi*2,100))])/2;
colormap(cmap')

% define time points for the movie frames
tidx = dsearchn(EEG.times',linspace(100,800,nframes)');
tidx = [tidx(1:end-1); tidx(end:-1:1)];

% azimuths for each frame
azs = round(linspace(0,360,length(tidx)))+90;

% create axes in specified locations
ax1_h = axes;
set(ax1_h,'Position',[.1 .1 .3 .8],'CameraViewAngle',6)

ax2_h = axes;
set(ax2_h,'Position',[.5 .1 .3 .8],'CameraViewAngle',6)


% octave: some features in headplotIndie are not (yet) implemented in
% Octave. You can try modifying the code to get it to work.

% create the movie!
for ti=1:length(tidx)
    
    axes(ax1_h);
    headplotIndie(squeeze(tf(:,1,tidx(ti))),'3dSplie.spl',[-3 3]);
    view(ax1_h,[-azs(ti) 20])
    
    axes(ax2_h);
    headplotIndie(squeeze(tf(:,2,tidx(ti))),'3dSplie.spl',[-3 3]);
    view(ax2_h,[azs(ti) 20])
    
    pause(.01)
end

%% end
