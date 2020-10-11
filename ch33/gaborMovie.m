

%% FMRI movie!!

%% load in fMRI data

filz = dir('*img');

for imgi=1:length(filz)
    
    tempdat = readnifti(filz(imgi).name);
    
    % initialize
    if imgi==1
        fmridat = zeros([ size(tempdat) length(filz) ]);
    end
    
    fmridat(:,:,:,imgi) = tempdat;
    
end

%%

clf

slice2plot = 35;
cdiscr = 64;

img2plot = squeeze(mean(fmridat(:,:,slice2plot,:),4))';
% normalize image to [0 64]
img2plot = img2plot - min(img2plot(:));
img2plot = cdiscr * img2plot./max(img2plot(:));


% zscore fmri brain map per voxel
fmridatZ = fmridat;
for ti=1:size(fmridatZ,4)
    tempvol = squeeze(fmridat(:,:,:,ti));
    fmridatZ(:,:,:,ti) = (tempvol-mean(tempvol(:))) ./ std(tempvol(:));
end
% fmridatZ = bsxfun(@minus,fmridat,mean(fmridat,4));
% fmridatZ = bsxfun(@rdivide,fmridatZ,std(fmridatZ,[],4));

% normalize image to [64 128]
fmridatZ = fmridatZ - min(fmridatZ(:));
fmridatZ = fmridatZ./max(fmridatZ(:));
fmridatZ = cdiscr+cdiscr*.5 + cdiscr*.5 * fmridatZ;

% same for stats map
stat2plot = squeeze(fmridatZ(:,:,slice2plot,1))';
stat2plot(stat2plot<110) = NaN;



h(1) = pcolor( img2plot );
hold on
h(2) = pcolor( stat2plot );
axis xy, axis image, axis off

colormap([ gray(cdiscr); hot(cdiscr) ]);
set(h,'linestyle','none')
set(gca,'clim',[0 128])

%

for ti=1:size(fmridatZ,4)
    
    stat2plot = squeeze(fmridatZ(:,:,slice2plot,ti))';
    stat2plot(stat2plot<107) = NaN;
    
    set(h(2),'Cdata',stat2plot);
    
    drawnow
    pause(.1)
end

%%




