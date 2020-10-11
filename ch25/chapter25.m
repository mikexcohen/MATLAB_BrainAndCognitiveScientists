%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 25
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% load and plot structural MRI data

% load structural mri file
[strMRI pixdim rotate dtype] = readnifti('MNI152_T1_1mm.nii');

% image difference slices
figure(1), clf
subplot(221)
imagesc(squeeze(strMRI(60,:,:))')
axis image, axis xy

subplot(222)
imagesc(squeeze(strMRI(:,180,:))')
axis image, axis xy

subplot(223)
imagesc(squeeze(strMRI(:,:,50))')
axis image, axis xy

colormap gray
set(gcf,'color','k')

% Note: If you look carefully, you can see the axes around the brains. 
% They are dark grey. You can turn them off with
axis off % this command applies only to the last-used axis

%% load in fMRI data

filz = dir('*img');

for imgi=1:length(filz)
    
    % read one volume
    tempdat = readnifti(filz(imgi).name);
    
    % initialize
    if imgi==1 % initialization only during first image!
        fmridat = zeros([ size(tempdat) length(filz) ]);
        volnums = zeros([ 1 length(filz) ]);
    end
    
    % put volume into larger matrix
    fmridat(:,:,:,imgi) = tempdat;
    
    
    % and get file name
    uscore = strfind(filz(imgi).name,'_');
    dotloc = strfind(filz(imgi).name,'.');
    volnums(imgi) = sscanf(filz(imgi).name(uscore+1:dotloc-1),'%g');
end

%% plot a few images and time courses

figure(2), clf, set(gcf,'color','w')

% image one slice (try some others!)
subplot(221)
imagesc(squeeze(fmridat(20,:,:,1))')
axis xy, axis image, axis off
set(gca,'clim',[0 12000])
title('Volume 20')

% same slice but later in time
subplot(222)
imagesc(squeeze(fmridat(20,:,:,end-10))')
axis xy, axis image, axis off
set(gca,'clim',[0 12000])
title([ 'Volume ' num2str(size(fmridat,4)-10) ])

colormap gray

% and a time course
subplot(212)
plot(squeeze(fmridat(20,40,40,:)))
title(sprintf('Time course of voxel I,J,K=%g,%g,%g',20,40,40))
xlabel('Time (volume)')
ylabel('Signal intensity (a.u.)')

%% experiment design

onsets = (6:12:84)+12;

timeline = zeros(length(filz),1);
for i=0:5
    timeline(onsets+i) = 1;
end

%% numerator of t-test

m0 = mean(fmridat(:,:,:,timeline==0),4);
m1 = mean(fmridat(:,:,:,timeline==1),4);

numerator = m1-m0;

%% denominator of t-test

v0 = var(fmridat(:,:,:,timeline==0),[],4);
v1 = var(fmridat(:,:,:,timeline==1),[],4);

denominator = sqrt( (v0/sum(timeline==0)) + (v1/sum(timeline==1)) );

%% t-statistic

% map of t-statistic values
tmap = numerator ./ denominator;

% fairly arbitrarily threshold at t<2.5
tthresh = tmap;
tthresh(tthresh<2.5) = 0;

% let's have a look at a thresholded slice
figure(3), clf
imagesc(squeeze(tthresh(:,:,35))')

%% threshold 1-voxel results

% slightly better threshold based on contiguously 
% significant regions. bwconncomp is in the image
% processing toolbox.
islands = bwconncomp(tthresh);
% find the number of pixels in each 
islandsizes = cellfun(@length,islands.PixelIdxList);

for ii=1:islands.NumObjects
    if islandsizes(ii)<3
        tthresh(islands.PixelIdxList{ii}) = 0;
    end
end

%% now plot

figure(4), clf

% pick a slice to show (fix to be axial)
slice2plot = 37;

% Color discritization. 64 is the Matlab standard.
% Try other numbers to see what happens to the plot. 
% Try, for example, 5 vs. 500
cdiscr = 64;

% extract slice and normalize image to [0 64]
img2plot = squeeze(m1(:,:,slice2plot))';
img2plot = img2plot - min(img2plot(:));
img2plot = cdiscr * img2plot./max(img2plot(:));


% same for stats map
stat2plot = squeeze(tthresh(:,:,slice2plot))';
stat2plot(stat2plot==0) = NaN;
% this time, normalize image to [64 128] because it will 
% be plotted on top of the structural map.
stat2plot = stat2plot - min(stat2plot(:));
stat2plot = cdiscr+cdiscr*.5 + cdiscr*.5 * stat2plot./max(stat2plot(:));

% uncomment the next line to reproduce the white voxels in the book figure
% stat2plot(stat2plot>1) = max(stat2plot(:));

% use pcolor function, which is similar to surf.
% Try replacing pcolor with surf and turning on rotate3d
h(1) = pcolor( img2plot );
hold on
h(2) = pcolor( stat2plot );
axis xy, axis image, axis off

colormap([ gray(cdiscr); hot(cdiscr) ]);
set(h,'linestyle','none')

%% end.
