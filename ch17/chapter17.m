%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 17
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% eigendecomposition

A = [3 1; 1 2];
[eigvecs,eigvals] = eig(A);

% define vectors
v1 = [.7 -.5]';
v2 = eigvecs(:,1);

% matrix-vector multiplication
v1A = A*v1;
v2A = A*v2;

% maximum value for plotting
xval = max([ abs(v1A); abs(v2A) ])*1.1;


figure(1), clf

subplot(131)
imagesc(A), axis square
title('Matrix A')


subplot(132)
plot([0 v1(1)],[0 v1(2)],'k','linew',4)
hold on
plot([0 v1A(1)],[0 v1A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Not an eigenvector!')



subplot(133)
plot([0 v2(1)],[0 v2(2)],'k','linew',4)
hold on
plot([0 v2A(1)],[0 v2A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Yes an eigenvector!')


%% PCA on simulated data

% data
x = [ 1*randn(1000,1) .4*randn(1000,1) ];

% rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data
y = x*R1;

% PCA of x (original data)
x = bsxfun(@minus,x,mean(x,1));
covmat = (x'*x) / (length(x)-1);
[evecsX,evalsX] = eig(covmat);

% PCA of y (correlated data)
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / (length(y)-1);
[evecsY,evalsY] = eig(covmat);


figure(2), clf
% plot original data
subplot(131)
plot(x(:,1),x(:,2),'m.','markersize',5)
set(gca,'xlim',[-5 5],'ylim',[-5 5])
hold on
plot(evalsX(1,1)*[0 evecsX(1,1)],evalsX(1,1)*[0 evecsX(2,1)],'k','linew',4)
plot(evalsX(2,2)*[0 evecsX(1,2)],evalsX(2,2)*[0 evecsX(2,2)],'k','linew',4)
xlabel('x-axis'), ylabel('y-axis')
axis square


subplot(132)
plot(y(:,1),y(:,2),'m.','markersize',5)
set(gca,'xlim',[-5 5],'ylim',[-5 5])
hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'k','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'k','linew',4)
xlabel('x-axis'), ylabel('y-axis')
axis square

% compute component scores
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

subplot(133)
plot(pc2,pc1,'m.')
set(gca,'xlim',[-5 5],'ylim',[-5 5])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square

%% vectors vs. values

x = rand(20);
covmat = (x'*x) / (length(x)-1);
[evecsX,evalsX] = eig(covmat); % note: mean-centering is not necessary here because this is not a covariance matrix

figure(3), clf
subplot(121)
imagesc(evecsX)
axis square

subplot(122)
imagesc(evalsX)
set(gca,'clim',[-.2 .2])
axis square

%%

figure(4), clf
hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'k','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'k','linew',4)
xlabel('x-axis'), ylabel('y-axis')
set(gca,'xlim',[-1 1],'ylim',[-1 1])


plot([0 evecsY(1,1)],[0 evecsY(2,1)],'r--','linew',2)
plot([0 evecsY(1,2)],[0 evecsY(2,2)],'r--','linew',2)
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')

set(gca,'xtick',-1:.25:1,'ytick',-1:.25:1)

%% eigenfaces! 

% read in one face to see what the data look like
dat = fread(fopen('faces/3500'));


figure(5), clf

% data are a vector... perhaps a plot?
subplot(211)
plot(dat)

% nope. image looks better
subplot(212)
imagesc(reshape(dat,128,128)')
axis image, axis off, colormap gray

%% import all faces

filz = dir('faces/*');
filz([filz.isdir]) = [];

allfaces = zeros(length(filz),length(dat));

for facei=1:length(filz)
    allfaces(facei,:) = fread(fopen([ 'faces/' filz(facei).name ]));
end

% reduce image size by masking
mask = false(128);
mask(40:100,20:120) = true;
maskdims = size(mask(any(mask,2),any(mask,1)));
allfaces = allfaces(:,mask);

%% PCA on faces

% mean-subtract
allfaces = bsxfun(@minus,allfaces,mean(allfaces,2));

% covariance
facecov = (allfaces'*allfaces) / length(allfaces);

% eigendecomposition (might take some seconds...)
[facevecs,facevals] = eig(facecov);

%% plot some eigfaces

figure(6), clf

for i=1:3
    subplot(2,3,i)
    imagesc(reshape(facevecs(:,end-i+1),maskdims(1),maskdims(2))')
    set(gca,'clim',[-.02 .02])
    axis image, axis off, title([ 'Eigenface component ' num2str(i) ])
end

for i=1:3
    subplot(2,3,i+3)
    imagesc(reshape(facevecs(:,i),maskdims(1),maskdims(2))')
    set(gca,'clim',[-.02 .02])
    axis image, axis off, title([ 'Eigenface component ' num2str(length(facevecs)+1-i) ])
end

% looks best in grayscale
colormap gray

%% reconstruct a face from PCs

% select a particular image to reconstruct
whichface = 10;

% specify how many components (i.e., dimensions)
numPCs = 10;

% compute the low-dimensional estimate
facescores = allfaces*facevecs(:,1:numPCs);


figure(7), clf

% plot the original (full-dimensional) face
subplot(121)
imagesc(reshape(allfaces(whichface,:),maskdims)')
axis square

% Here we want to plot the low-dimensional reconstruction 
% of the face from the PCs, but somehow it looks awful.
% Check my code -- did I make a mistake? (Ok, ok, obviously
% there is an error somewhere. It's your job to find and fix it!)
subplot(122)
imagesc(reshape(facescores(whichface,:)*facevecs(:,1:numPCs)',maskdims)')
axis square


colormap gray

%% ICA vs. PCA

% generate data

% data
x = [ 1*randn(1000,1) .05*randn(1000,1) ];

% rotation matrix
th = -pi/6;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];
th = -pi/3;
R2 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data
y = [ x*R1 ; x*R2 ];


figure(8), clf
subplot(221)
plot(y(:,1),y(:,2),'o')

datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')


% now PCA
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / length(y);
[evecsY,evalsY] = eig(covmat);

hold on
plot([0 evecsY(1,1)],[0 evecsY(2,1)],'color',colorrgb(4,:),'linew',4)
plot([0 evecsY(1,2)],[0 evecsY(2,2)],'color',colorrgb(5,:),'linew',4)


subplot(222)
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

plot(pc2,pc1,'ms')
datarange = max([pc1(:); pc2(:)])*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square
title('Data in PC space')





% now ICA
subplot(223)
plot(y(:,1),y(:,2),'o')
datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])

ivecs = jader(y');
hold on
plot([0 ivecs(1,1)],[0 ivecs(2,1)],'color',colorrgb(4,:),'linew',4)
plot([0 ivecs(1,2)],[0 ivecs(2,2)],'color',colorrgb(5,:),'linew',4)
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')



subplot(224)
ic_scores = ivecs*y';
plot(ic_scores(1,:),ic_scores(2,:),'ms')
datarange = max(ic_scores(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('IC1 axis'), ylabel('IC2 axis')
axis square
title('Data in IC space')

%% end.
