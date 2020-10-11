%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 12
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% "manual" time-domain convolution

% signal
signal = zeros(1,30);
signal(11:19) = 1;

% kernel (Gaussian)
kernel = [0 .31 .77 1 .77 .31 0];

% normalize the kernel (try uncommenting to see 
% the effect on the result)
kernel = kernel./sum(kernel);


% sizes
N = length(signal);
M = length(kernel);
halfKern = floor(M/2);

% Data that we'll use for convolution (must be zero-padded).
dat4conv = [ zeros(1,M-1) signal zeros(1,M-1) ];

% initialize convolution output
conv_res = zeros(1,N+M-1);

% run convolution
for ti=M+1:N-M+1
    tempdata = dat4conv(ti:ti+M-1);
    
    % compute dot-product (don't forget to flip the kernel backwards!)
    conv_res(ti) = sum( tempdata.*kernel(end:-1:1) );
end

% trim the wings
conv_res = conv_res(halfKern+1:end-halfKern);

figure(1), clf

% plot signal
subplot(311)
plot(signal)
set(gca,'xlim',[0 N],'ylim',[min(signal)-.05 max(signal)+.05])

% plot kernel
subplot(312)
plot(kernel)
set(gca,'xlim',[0 N])

% plot convolution result
subplot(313)
plot(conv_res,'o-','linew',2,'markerface','g','markersize',9)
set(gca,'xlim',[0 N],'ylim',[min(signal)-.05 max(signal)+.05])

% Run this next line to show that our manual convolution 
% matches the output of Matlab's conv function.
%hold on, plot(conv(signal,kernel,'same'),'r*')

%% time-domain convolution as frequency-domain multiplication

% step 1: get the right sizes
nData = length(signal);
nKern = length(kernel);
nConv = nData+nKern-1;
nHfkn = floor(nKern/2);

% step 2: FFTs
signalX = fft(signal,nConv);
kernelX = fft(kernel,nConv);
hold on, plot(conv(signal,kernel,'same'),'r*')

% step 3: normalize the kernel
kernelX = kernelX./max(kernelX);

% step 4: point-wise multiplication
convres = ifft( signalX.*kernelX );

% could be broken down to 4a and 4b:
%convres = signalX.*kernelX;
%convres = ifft( convres );

% step 5: trim wings
convres = convres(nHfkn+1:end-nHfkn);

hold on
plot(convres,'rd','linew',4)

%% convolution between signal and noise

% signal of random noise
signal = randn(1,1000);

% gaussian kernel (width is hard-coded)
x = -2:.01:2;
kernel = exp( (-x.^2) / .1 );


% step 1: get the right sizes
nData = length(signal);
nKern = length(kernel);
nConv = nData+nKern-1;
nHfkn = floor(nKern/2);

% step 2: FFTs
signalX = fft(signal,nConv);
kernelX = fft(kernel,nConv);

% step 3: normalize the kernel
kernelX = kernelX./max(kernelX);

% step 4: point-wise multiplication
convres = ifft( signalX.*kernelX );

% could be broken down to 4a and 4b:
%convres = signalX.*kernelX;
%convres = ifft( convres );

% step 5: trim wings
convres = convres(nHfkn+1:end-nHfkn);

figure(6), clf
subplot(231)
plot(signal)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: signal')

subplot(232)
plot(kernel)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: kernel')

subplot(233)
plot(convres)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: convolution result')


subplot(234)
plot(2*abs(signalX))
set(gca,'xlim',[0 nConv])
xlabel('Frequency (a.u.)')
title('Frequency domain: signal')

subplot(235)
plot(2*abs(kernelX))
set(gca,'xlim',[0 80])
xlabel('Frequency (a.u.)')
title('Frequency domain: kernel')

subplot(236)
plot(2*abs(signalX.*kernelX))
set(gca,'xlim',[0 80])
xlabel('Frequency (a.u.)')
title('Frequency domain: convolution result')

%% 2D convolution

% import picture (We'll use only the red layer, because today feels like a red day, don't you think?)
pic = imread('saturn.png');
pic = double(squeeze(pic(:,:,1)));

% create gaussian
[x,y]  = meshgrid(-250:250);
s = 30;
gaus2d = exp( -(x.^2 + y.^2)/(2*s^2) );


% step 1
N     = size(pic);
M     = size(gaus2d);
nConv = N+M-1;
halfK = floor(M/2);

% step 2
picX  = fft2(pic,nConv(1),nConv(2));
gausX = fft2(gaus2d,nConv(1),nConv(2));

% step 3
gausX = gausX ./ max(gausX(:));

% step 4
cr = ifft2( picX.*gausX );

% step 5
cr = cr( halfK(1)+1:end-halfK(1) , halfK(2)+1:end-halfK(2) );


figure(2), clf
subplot(231)
imagesc(pic)
axis square, axis off
title('Saturn in space')

subplot(232)
imagesc(gaus2d)
axis square, axis off
title('Gaussian in space')

subplot(233)
imagesc(cr)
axis square, axis off
title('Convolution result in space')


subplot(234)
imagesc(fftshift(log10(abs(picX))))
set(gca,'clim',[2.5 5])
axis square, axis off
title('Saturn in frequency')

subplot(235)
imagesc(fftshift(log10(abs(gausX))))
set(gca,'clim',[-20 0])
axis square, axis off
title('Gaussian in frequency')

subplot(236)
imagesc(fftshift(log10(abs( picX.*gausX ))))
set(gca,'clim',[-20 0])
axis square, axis off
title('Convolution result in frequency')

%% end.
