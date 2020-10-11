%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 27
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% load in the Saturn picture

saturn = double( imread('saturn.png') )./255;
% imagesc(saturn)

% first, add some noise
noiselevel = .5;
saturn = saturn + randn(size(saturn))*noiselevel;

% let's see how noisy it is
figure(1), clf
imagesc(saturn), axis square

%% mean-smoothing the wrong way

% What are the problems with this code?

k=9;
for i=k+1:size(saturn,1)
    for j=k+1:size(saturn,2)
        temp = saturn(i-k:i+k,j-k:j+k);
        saturn(i,j) = mean(temp(:));
    end
end

%% mean smoothing the right way

k=9;

saturn = double( imread('saturn.png') )./255;
saturnFilt = saturn;

for dimi=1%:size(saturn,3)
    for i=k+1:size(saturn,1)-k
        for j=k+1:size(saturn,2)-k
            temp = saturn(i-k:i+k,j-k:j+k,dimi);
            saturnFilt(i,j,dimi) = mean(temp(:));
        end % end loop over rows
    end % end loop over columns
end % end loop of dimensions

figure(2), clf
subplot(121)
imagesc(squeeze(saturn(:,:,1)))
axis xy, axis square, axis off
set(gca,'clim',[-1 1])
title('Original')


subplot(122)
imagesc(squeeze(saturnFilt(:,:,1)))
axis xy, axis square, axis off
set(gca,'clim',[-1 1])
title([ 'Mean-smoothed with k=' num2str(k) ])


%% add large-amplitude spikes

% we need to load in the image again because I didn't save the original
% before adding Gaussian noise to it (bad programming!).
saturn = double( imread('saturn.png') )./255;

% let's just use the first color dimension for simplicity
saturn = squeeze(saturn(:,:,1));

% add noise to a random subset of pixels
spikelocs = randsample(1:numel(saturn),round(.1*numel(saturn))); % for octave: randsample in the statistics package
noisySaturn = saturn;
noisySaturn(spikelocs) = 123456789;


figure(3), clf, colormap bone

% show the original
subplot(131)
imagesc(saturn)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title('Original')

% show the noisified version
subplot(132)
imagesc(noisySaturn)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title('Noisy image')



% k parameter for median filter
k=9;

% notice how I use different variables to avoid overwriting
saturnFilt = noisySaturn;

for i=k+1:size(saturn,1)-k
    for j=k+1:size(saturn,2)-k
        temp = saturn(i-k:i+k,j-k:j+k);
        saturnFilt(i,j) = median(temp(:));
    end % end loop over rows
end % end loop over columns


% and plot the filtered result. Looks pretty good!
subplot(133)
imagesc(saturnFilt)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title([ 'Median-smoothed with k=' num2str(k) ])

%% Gaussian-based image smoothing

% parameters for Gaussian
gx = -20:20;
gaus2d = zeros(length(gx));

% equal width in both directions. We might call this an isotropic Gaussian.
sx = 5;
sy = 5;

% create the Gaussian point-by-point
for xi=1:length(gx)
    for yi=1:length(gx)
        gaus2d(xi,yi) = exp( -(  (gx(xi)^2)/(2*sx^2)  +  (gx(yi)^2)/(2*sy^2) ));
    end
end

% plot the 2D gaussian
figure(4), clf
imagesc(gaus2d)
axis image
% hint: try using surf. It looks a bit like a... traffic cone?



% Now smooth the spike-field coherence image.
% (Hint: copy code from previous chapters.)
figure(5), clf
subplot(121)
imagesc(timevec,[],spikeLFP)
set(gca,'clim',[-1000 1000])

subplot(122)
smo = conv2(spikeLFP,gaus2d,'same');
imagesc(timevec,[],smo./sum(gaus2d(:)))
set(gca,'clim',[-500 500])

%% filter Saturn 

% load in a fresh copy, just in case
saturn = double( imread('saturn.png') )./255;


% get FFT (remember: 2D!). For simplicity,
% we use only the first color dimension.
saturnX = fft2(squeeze( saturn(:,:,1) ));

% get sizes of image and midpoints
imgdims = size(saturnX);
midX = round(imgdims(2)/2);
midY = round(imgdims(1)/2);

% size of the filter in pixels.
nPix2use = 100;

% create low-pass filter kernel
loPass2d = zeros(imgdims(1:2));
loPass2d(midY-nPix2use:midY+nPix2use,midX-nPix2use:midX+nPix2use) = 1;

% create high-pass filter kernel
hiPass2d = ones(imgdims(1:2));
hiPass2d(midY-nPix2use:midY+nPix2use,midX-nPix2use:midX+nPix2use) = 0;


% and let's see how it looks. First, the (shifted) power spectrum
figure(6), clf
subplot(221)
imagesc(fftshift(log(abs(saturnX))));
axis off, axis square
title('Full power spectrum')


% next, the filter
subplot(222)
imagesc(loPass2d)
axis off, axis square
title('Low-pass filter')

% apply the low-pass filter
filtimg = real(ifft2( saturnX.*fftshift(loPass2d) ));
subplot(223)
imagesc(filtimg)
axis off, axis square
title('Low-pass image')

% and the high-pass filter
filtimg = real(ifft2( saturnX.*fftshift(hiPass2d) ));
subplot(224)
imagesc(filtimg)
axis off, axis square
set(gca,'clim',[-.2 .2]/10)
title('High-pass image')

%% end
