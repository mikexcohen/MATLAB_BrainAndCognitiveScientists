%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 26
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% value-based thresholding

% load time-frequency data
load tfmat.mat

% before running the code, check out the contents of this file
whos

% define binarized threshold map
threshmap = p<.05;

% color limit for the image
clim = [-3 3];

% plot the raw image
figure(1), clf
subplot(221)
contourf(tf,40,'linecolor','none')
set(gca,'clim',clim)
xlabel('Time (a.u.)'), ylabel('Frequency (a.u.)')


% plot the thresholded and binarized map
subplot(222)
contourf(threshmap,1,'linecolor','none')
xlabel('Time (a.u.)'), ylabel('Frequency (a.u.)')


% show contours on the raw image
subplot(223)
contourf(tf,40,'linecolor','none')
hold on
contour(threshmap,1,'linecolor','k','linew',2)
set(gca,'clim',clim)
xlabel('Time (a.u.)'), ylabel('Frequency (a.u.)')


% segment the image
tfthresh = tf;
tfthresh(~threshmap) = 0;

subplot(224)
contourf(tfthresh,40,'linecolor','none')
set(gca,'clim',clim)
xlabel('Time (a.u.)'), ylabel('Frequency (a.u.)')

%% illusration of bwconncomp output

% Inspect the fields of this output:
% What information is contained in each field?
islands = bwconncomp(threshmap);

% initialize matrices
stepmap  = zeros(islands.ImageSize);
avepower = zeros(islands.NumObjects,1);

for i=1:islands.NumObjects
    % populate the map with numbered blobs for visualization
    stepmap(islands.PixelIdxList{i}) = i;
    
    % get the average power from the TF map for each blob (island)
    avepower(i) = mean(tf(islands.PixelIdxList{i}));
end

figure(2), clf
imagesc(stepmap)
axis xy

%% import the MRI image

% this function is not standard, and is located in a different folder.
% The easiest way to run it is to add that folder to the path.
addpath('../ch25')
smri = readnifti('../ch25/MNI152_T1_1mm.nii');
% Question: Do you need to specify the location of the file? Why or why not?


figure(3), clf
subplot(311)
% WARNING: Depending on your computer's graphics card, 
% this next line might take a while to render.
% hist(smri)
xlabel('Value'), ylabel('Count')

subplot(312)
hist(nonzeros(smri),10000)
xlabel('Value'), ylabel('Count')

subplot(313)
hist(log10(nonzeros(smri)),10000)
xlabel('Value'), ylabel('Count')

%% segmenting the MRI

% based on visual inspection of the log10 histogram
threshold1 = 10.^[3.0 3.5];
threshold2 = 10.^[3.7 3.8];

smriThresh1 = smri>threshold1(1) & smri<threshold1(2);
smriThresh2 = smri>threshold2(1) & smri<threshold2(2);


% show two slices that are at the lower threshold
figure(4), clf
subplot(221)
imagesc(squeeze(smriThresh1(60,:,:))')
axis image, axis xy, axis off

subplot(222)
imagesc(squeeze(smriThresh1(:,180,:))')
axis image, axis xy, axis off



% show the sane two slices at the higher threshold
subplot(223)
imagesc(squeeze(smriThresh2(60,:,:))')
axis image, axis xy, axis off

subplot(224)
imagesc(squeeze(smriThresh2(:,180,:))')
axis image, axis xy, axis off



colormap gray
set(gcf,'color','k')

%% the next image to segment

% calcium image take from sample dataset:
% http://dylan-muir.com/projects/focusstack_stimserver/
load image4imageseg

figure(5), clf
imagesc(im)

colormap gray
axis image
colorbar

%% segment the image

% 120 is an arbitrary threshold from trial-and-error guessing
% note that here we use bwlabel instead of bwconncomp. What's the
% difference?
[islands,numblobs] = bwlabel(im>120);

% reconstruct the bwconncomp output.
pixelIdxList = cell(1,numblobs);
for i=1:numblobs
    pixelIdxList{i} = find(islands==i);
end


figure(6), clf

subplot(221)
imagesc(islands) % what happened to the blobs on the left of the image?
colormap gray
axis image

subplot(222)
imagesc(log(islands)) % does the log of the values help?
axis image

subplot(223)
imagesc(logical(islands))
axis image

% get rid of blobs with fewer than 21 continuous pixels.
% Also an arbitrary threshold from trial-and-error guessing.
% Try changing the threshold to see the effect.
for i=1:numblobs
    if sum(islands(:)==i)<21
        islands(islands==i)=0;
    end
end

% and the new thresholded image
subplot(224)
imagesc(logical(islands))
axis image

%% get the average luminance value from all clusters

activity = zeros(1,numblobs);

for i=1:numblobs
    activity(i) = mean(im(islands==i));
end

figure(7), clf
bar(activity)
xlabel('Island number'), ylabel('Mean activity level')

% but how many clusters are there really?

%% grid discretization

% specify number of points and number of discritizations
n = 300;
k =  7;

% convert indices to grids
[gridX,gridY] = ndgrid(1:n,1:n);

% let's see how these outputs look
figure(8), clf
subplot(231), imagesc(gridX), colorbar, axis image
subplot(232), imagesc(gridY), colorbar, axis image

% discretization
gridX = ceil( k*gridX./n );
gridY = ceil( k*gridY./n );

% and let's see how these look
subplot(234), imagesc(gridX), colorbar, axis image
subplot(235), imagesc(gridY), colorbar, axis image

% scale up, combine, scale down (why?)
tempG = gridX + gridY*(k*1000);
u=unique(tempG);
grids = zeros(n);

for ui=1:length(u)
    grids(tempG==u(ui)) = ui;
end

subplot(236), imagesc(grids'), colorbar, axis image

%% Sierpinksi triangle, in dots and in image

% parameters...
N = 10000; % number of points in the image
ds_factor = 40; % down-sampling factor for analysis
Nds = round(N/ds_factor); % convert the down-sampling factor

% further initializations...
[sx,sy] = deal( zeros(1,N) );
siertri = zeros(Nds);

% populate the Sierpinksi image
for i = 2:N
    k = ceil(rand*3);
    sx(i) = sx(i-1)/2 + (k-1)*.25;
    sy(i) = sy(i-1)/2 + (k==2)*.5;
end

% here is the image as a bunch of points
figure(9), clf
plot(sx,sy,'m^')
axis off
set(gcf,'color','k')


% and the image as an image
figure(10), clf
wherezeros = sx==0 | sy==0;
sx(wherezeros)=[]; sy(wherezeros)=[];
siertri(sub2ind([Nds,Nds],ceil(Nds*sy),ceil(Nds*sx))) = 1;
imagesc(siertri), colormap gray, axis xy

%% box-counting analysis of sierpinksi image

% define a "mother" discrization grid (O for original)
[OgridX,OgridY] = ndgrid(1:Nds,1:Nds);

% discritizations (box sizes). Why does it count backwards? Does that matter?
ks2use = 50:-2:4;
tot = zeros(size(ks2use));

% loop over image sizes
for ki=1:length(ks2use)
    
    k = ks2use(ki);
    
    % downsample 'mother' grid lines
    
    % further discretization and differentiation
    % to get unique values for each box
    
    
    % find whether each box contains part of the image
    for ui=1:length(u)
        
    end
    
    % sum over all boxes from the previous loop
    tot(ki) = ;
end

% then plot the results
figure(12), clf
imagesc(2-siertri - 1*(tempG==u(14)))
axis xy, colormap gray

figure(11), clf
plot(log10(ks2use),log10(tot),'o-')
xlabel('log10(box size)')
ylabel('log10(image size)')
set(gca,'xlim',[.5 1.8])

%% end.
