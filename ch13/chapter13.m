%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 13
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% intro to interpolation

% 'measured' data
data      = [1 4 3 6 2 19];
datatimes = 1:6;

% requested (interpolated) time points
newtimes  = 1:.001:6;

% different interpolation options
interpOptions = {'linear';'spline';'next';'nearest'};

figure(1), clf
for methodi=1:length(interpOptions)
    
    % define interpolation object
    % this function is not implemented in Octave; see interp1)
    F = griddedInterpolant(datatimes,data,interpOptions{methodi});
    
    % query that object at requested time points
    newdata = F(newtimes);
    
    % and plot
    subplot(2,2,methodi)
    plot(newtimes,newdata,'k','linew',3)
    hold on
    plot(datatimes,data,'go','markersize',15,'markerfacecolor','m')
    set(gca,'xlim',[.5 6.5])
    title([ '''' interpOptions{methodi} '''' ]) % 4 single quotes here to get a single quote of text!
end

%% interpolation vs. extrapolation

% which row to interpolate?
% Try 7 (interpolation) or 17 (extrapolation)
point2interp = 17;

% get the landscape
[A,B,C] = deal( peaks(10) );

% define a regular grid but missing row 7. We will then interpolate
% or extrapolate using linear and spline methods.
[x,y] = ndgrid([1:6 8:10],1:10);
Flin  = griddedInterpolant(x,y,A([1:6 8:10],:),'linear');
Fspl  = griddedInterpolant(x,y,A([1:6 8:10],:),'spline');

% query data points
B(7,:) = Flin(repmat(point2interp,1,10),1:10);
C(7,:) = Fspl(repmat(point2interp,1,10),1:10);

figure(2), clf
subplot(221)
imagesc(A)
set(gca,'clim',[-7 7])
title('Original')

subplot(222)
imagesc(B)
set(gca,'clim',[-7 7])
title('linear')

subplot(223)
imagesc(C)
set(gca,'clim',[-7 7])
title('spline')

% plot 
subplot(224)
plot(1:10,A(7,:),'k'), hold on
plot(1:10,B(7,:),'b')
plot(1:10,C(7,:),'r')
set(gca,'xlim',[.8 10.2])

legend({'original';'linear';'spline'})

%% an aside on ndgrid
% you'll learn more about making grids and using them for image
% processing in Chapter 26. For now, try to get a feel for what
% the outputs are by running each line and inspecting the outputs.

% no difference...
x = ndgrid(1:3);
x = ndgrid(1:3,2);

% also no difference
x = ndgrid(1:3,[1 2]);
x = ndgrid(1:3,[1 3]);



% now with 2 outputs
[x,y] = ndgrid(1:3);
[x,y] = ndgrid(1:3,2);

% What's the difference between these lines?
[x,y] = ndgrid(1:3,[1 2]);
[x,y] = ndgrid(1:3,[1 3]);


figure(3), clf
[x,y] = ndgrid(1:30,[1:10 21:30]);
subplot(221)
plot(1:20,x), axis square

subplot(223)
imagesc(x), axis square

subplot(222)
plot(y,1:30), axis square

subplot(224)
imagesc(y), axis square

%% 2D interpolation using real EEG data

% load some sample EEG data
load EEGexample.mat

% convert polar to cartesian coordinates
[eX,eY] = pol2cart(pi/180*[chanlocs.theta],[chanlocs.radius]);

% interpolation factor, and define grid spacing
intFact = 100;
interpX = linspace(min(eX),max(eX),intFact);
interpY = linspace(min(eY),max(eY),intFact);

% now define grid in which to interpolate
[gridX,gridY] = ndgrid(interpX,interpY);


% shall we have a look at the grids?
figure(4), clf
subplot(222)
imagesc(gridX)

subplot(223)
imagesc(gridY)

subplot(224)
imagesc(gridX+gridY)


% now define the interpolation object, this time using scatteredInterpolant
F = scatteredInterpolant(eX',eY',eeg(:,300),'linear','none');
interpDat = F(gridX,gridY);


figure(5), clf

% image the interpolated map
subplot(121)
imagesc(interpX,interpY,interpDat')
axis off, axis image
axis([-.7 .6 -.6 .6])

% draw actual electrode positions
hold on
plot(eX,eY,'ko','markerfacecolor','w','markersize',8)


% this is what the 'map' would look like without interpolation
subplot(122)
scatter(eX,eY,50,eeg(:,300),'filled')
axis image
set(gca,'color','k')
axis([-.7 .6 -.6 .6])

%% using interp1

% Define the 'real' data. Similar to the first cell.
data = [1 4 3 6 2 19];
datatimes = 1:6;
newtimes = 1:.001:6;

% test four options (linear and spline work in Octave)
interpOptions = {'linear';'spline';'next';'nearest'};

figure(6), clf
for methodi=1:length(interpOptions)
    
    % get the new data
    newdata = interp1(datatimes,data,newtimes,interpOptions{methodi});
    
    % specify subplot
    subplot(2,2,methodi)
    
    % plot data
    plot(newtimes,newdata,'k','linew',3)
    hold on
    plot(datatimes,data,'go','markersize',15,'markerfacecolor','m')
    
    % adjustments..
    set(gca,'xlim',[.5 6.5])
    legend({'interpolated';'original'})
end

%% zero-padding

origN =  10;
padN  = 100;

% random data
data = randn(origN,1);

% padded data
datapad = ifft( fft(data)/origN ,padN)*padN;


% and plot!
figure(7), clf
plot(linspace(0,1-1/length(datapad),padN),real(datapad),'bo-','markersize',10)
hold on
plot(linspace(0,1-1/origN,origN),data,'rs-','markersize',10)

legend({'zero-padded';'original'})

%% downsample

origSrate = 1000;
newSrate  =  250;

% specify time vector
time = 0:1/origSrate:1-1/origSrate;

% smooth the data with a Gaussian
data = conv( randn(length(time),1), gausswin(400,20) ,'same');

dsFact = origSrate/newSrate;


datads = data(1:dsFact:end);
timeds = time(1:dsFact:end);


figure(8), clf
plot(time,data,'b.-','markersize',10)
hold on
plot(timeds,datads,'ro')

legend({'original';'downsampled'})

%% end
