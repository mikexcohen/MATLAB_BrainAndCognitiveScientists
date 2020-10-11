%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 9
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% creating and destroying figures

% open 3 figures
figure, figure, figure


% specify the figure number to open
figure(1) 
figure(10)
fignum = 400;
figure(fignum)


close % closes only active figure
close(10) % closes figure #10
close([2 3 5:8])
close all

%% basic line plotting

% very basic (though surprising useful)
plot(1:10)

% slightly less basic
x = 1:2:20;
y = x.^2;
plot(x,y)

% demonstrating the 'hold' functionality
plot(x,y)
hold on
plot(x,log(y))
hold off
plot(x,y.^(1/3))


cla % cla = clear axis
plot(x,y/50,'r')
hold on
plot(x,log(y),'k')
plot(x,y.^(1/3),'m')

hold off % now the second line will overwrite the following line
plot(x,y,'ro-','linewidth',9) % default width is 1
plot(x,y,'ro-','linewidth',9,'markersize',100) % bigger markers (too big?)

%% legend

% I like to control figure numbering.
figure(1), clf

plot(x,y/40,'bp')
hold on
plot(x,log(y),'r*--')

% how does the order of the legend items 
% relate to the order of the plotting?
legend( {'y=x^2/40' ; 'y=log(x^2)'} )

%% bar plots and errorbar plots

figure(2), clf

% basic bar plot
bar(x,y)

% specify the width of the bars
bar(x,y,.2) % try other numbers


% error bars
e = 100*rand(size(x));
errorbar(x,y,e) % symmetric bars
errorbar(x,y,e/2,e/8) % asymmetric bars

% combination of bars and errorbars
bar(x,y)
hold on
errorbar(x,y,e,'.')

%% scatterplot

% basic scatterplots are the same as the 'plot' function
figure(3), clf
scatter(x,y,'o')

% with more options, the scatterplot function can be used
% to show multidimensional results

% create some data
n      = 100;
frate  = linspace(10,40,n) + 10*rand(1,n); % firing rates
fvar   = frate + 5*randn(1,n);
ndepth = linspace(100,1000,n);

% color can be used as a third dimension of information
scatter(frate,fvar,100,ndepth,'filled') % Octave: try 10 instead of 100!
% what does the '100' mean in the scatter function input?

%% histograms

nRand = 10000; % soft-coding!

figure(4), clf

% basic histogram function
r = randn(nRand,1);
hist(r,50) % 50 bins (default is 10)


% slightly more advanced use of hist function
ru = rand(nRand,1);
[y_r,x_r] = hist(r,50); % y outputs first
[y_ru,x_ru] = hist(ru,50);

plot(x_r,y_r,'k','linew',2) % linew = linewidth
hold on
plot(x_ru,y_ru,'r','linew',2)
legend({ 'randn';'rand' })
xlabel('Value'), ylabel('Count') % do these two functions need to be on the same line?

%% subplots

figure(5), clf

subplot(221) % commas unnecessary for <10 subplots
plot(r,'k')
title('normal random numbers')
xlabel('Indices'), ylabel('Values')

subplot(2,2,2) % but commas can improve readability
plot(ru,'k')
title('uniform random numbers')
xlabel('Indices'), ylabel('Values')

subploti = 3;
subplot(2,2,subploti) % commas are necessary when using variables to index a subplot
plot(x_r,y_r,'k')
title('distribution of normal')
xlabel('Values'), ylabel('Count')

subplot(224)
plot(x_ru,y_ru,'k')
title('distribution of uniform')
xlabel('Values'), ylabel('Count')


% mixing geometry
clf
subplot(221), plot(r) % terse code can be stacked
subplot(222), plot(ru)
subplot(212), plot(x_r,y_r) % note the difference in subplot input
hold on, plot(x_ru,y_ru,'r')

%% going crazy with subplots in one figure!

figure(6), clf

% The functions 'text' and 'set' might be new to you.
% For now, just pay attention to the subplot geometry.
txt_h = zeros(1,4);

subplot(311),     txt_h(1) = text(.5,.5,'subplot(311)');
subplot(345),     txt_h(2) = text(.5,.5,'subplot(345)');
subplot(246),     txt_h(3) = text(.5,.5,'subplot(246)');
subplot(739),     txt_h(4) = text(.5,.5,'subplot(739)');
subplot(224),     txt_h(5) = text(.5,.5,'subplot(224)');
subplot(8,8,57),  txt_h(6) = text(.5,.5,'subplot(8,8,57)');


% center align all text objects
set(txt_h,'HorizontalAlignment','Center')

%% patches

% define two vectors
x = [1 2 3 4 3 2 1];
y = [9 9 7 4 1 3 2];

figure(7), clf

% plot edge markers
plot(x,y,'o','markerfacecolor','g','markersize',15)
h=patch(x,y,'r')


% note that patches are held on by default!
[~,idx] = sort(x);
patch(x(idx),y(idx),'r')


% Let's try it again in a new figure.
% Now you see that changing the order of the
% coordinates changes the patch.
figure(8), clf
patch(x(idx),y(idx),'r')
hold on
plot(x,y,'o','markerfacecolor','g','markersize',15)

%% images

figure(9), clf

% read in saturn picture (comes with Matlab)
pic = imread('saturn.png');
imagesc(pic)

% try different axis properties:
axis image
% axis square
% axis normal


% oops (but why?)
pic2 = pic;
pic2(:,:,4) = pic(:,:,1);
imagesc(pic2)


% showing the different coError using image
colorchans = { 'red';'green';'blue' };
for chani=1:3
    subplot(2,2,chani+1)
    imagesc(pic(:,:,chani))
    axis off
    set(gca,'clim',[0 255])
    title([ colorchans{chani} ' channel' ])
end

%% contour plots

% if you have a slow graphics card or are running Octave, 
% be patient with this section...

figure(10), clf
contourf(pic(:,:,1))


% trying different settings of contour (run each line separately)
contourf(pic(:,:,1),40) % 10 is the default number 
contourf(pic(:,:,1),40,'linecolor','m')
contourf(pic(:,:,1),40,'linecolor','none') % my favorite option

% repeat using the contour function (no fillings)
contour(pic(:,:,1))
contour(pic(:,:,1),1,'linecolor','k'), axis off
contour(pic(:,:,1),'linecolor','m'), axis off
set(gcf,'color','k') % you'll learn about set soon

%% surf

figure(11), clf
surf(pic(:,:,1))

% now a bit better
shading interp
axis off

%% contourf vs. imagesc

figure(12), clf
subplot(121)
contourf(pic(:,:,1),40,'linecolor','none') 
title('contourf')

subplot(122)
imagesc(pic(:,:,1))
title('imagesc')

% what are the differences between these plots?

%% get

figure(13), clf
plot(rand(3))

% basic get usage
get(gca,'xlim')
yTik = get(gca,'ytick');

% list all properties
get(gca)

% using get to inform plotting
hold on
plot( get(gca,'xlim') ,[.1 .6],'k--','linew',3)

%% set

% 'get' accesses; 'set' changes

set(gca,'xlim',[0 2])
set(gca,'ytick',[0 .5 .8 .91])

% many changes in one set function
% note the pattern of value-property
set(gca,'xlim',[0 2],'ytick',[.5 1],'xtick',0:.25:2)

%% beyond gca: using plot handles

close all % sometimes it's good to clear the mental workspace

figure(14), clf

% the output of plot (and many other plotting functions) 
%   is a handle to that plot object
line_h = plot(1:10,(1:10).^2);

get(line_h)
get(gca) % different from the previous line!

set(line_h,'linewidth',4,'marker','o','markeredgecolor','k')

% now close the figure
close
% and try again:
set(line_h,'linewidth',4,'marker','o','markeredgecolor','k')
% uh oh... (why?)

%% many plot objects, many handles

% initialize handles vector
plot_hs = zeros(1,100);

% draw lots of lines with different handles
figure(15), clf
hold on
for i=1:100
    plot_hs(i) = plot( randn( max(1,round(rand*10)) ,1) );
end

% now you can access each or many plot objects
set(plot_hs(1:50),'color','k')
set(plot_hs(25:75),'marker','o')
set(plot_hs([1:10 20:5:100]),'linewi',4)

%% gcf instead of gca

set(gcf,'color','m','name','Results of experiment 2b')

%% text

figure(16), clf
text(.6,.4,'Yo!')

% just in case you didn't get the memo
for i=1:1000, text(rand,rand,'Yo!'), end


% plot area is not adjusted for new text objects
text(1.05,.7,'outside')


% using set to change text properties
clf
txt_h = text(.5,.5,'Hello')

set(txt_h,'Position',[.2 .7])
set(txt_h,'color','m','String','Zoidberg')
set(txt_h,'HorizontalAlignment','Center','FontSize',30)

% Special characters are used for super/subscript and Greek letters
set(txt_h,'String','Z_oidbe^r^g is an \alpha\_crab')
ylabel('Power (\muV^2)')

%% colormaps

jetmap = jet;
bonmap = bone;

% define my own colormap
circmap = (1+[cos(linspace(0,pi*2,100)); sin(linspace(0,pi*2,100)); cos(linspace(0,pi*2,100) + pi)])/2;


figure(17), clf
subplot(311)
plot(jetmap,'linew',2)
set(gca,'ylim',[-.1 1.1],'xlim',[0 length(jetmap)+1])
legend({'R';'G';'B'})
xlabel('image value range'), ylabel('Color intensity')
title('JET')

subplot(312)
plot(bonmap,'linew',2)
set(gca,'ylim',[-.1 1.1],'xlim',[0 length(bonmap)+1])
legend({'R';'G';'B'})
xlabel('image value range'), ylabel('Color intensity')
title('BONE')

subplot(313)
plot(circmap','linew',2)
set(gca,'ylim',[-.1 1.1],'xlim',[0 length(bonmap)+1])
legend({'R';'G';'B'})
xlabel('image value range'), ylabel('Color intensity')
title('CIRCULAR')

%% let's see colormaps in action

% let's use a 2D Gaussian. 
lims = [-31 31];
[x,y] = ndgrid(lims(1):lims(2),lims(1):lims(2));
gaus2d = exp( -(x.^2 + y.^2) / 200 );

% image the gaussian
figure(18), clf
imagesc(gaus2d)


% try different colormaps
colormap jet % the default in neuroscience
colormap bone
colormap summer
colormap hsv

% you can also define the resolution of the colormap
%  (doesn't work in octave...)
colormap jet(3)
colormap jet(300)
colormap jet(64) % 64 is the default


% now let's make our own simple colormap
myCmap = repmat( linspace(0,1,20)' ,1,3);
colormap(myCmap)
% Does inspecting myCmap tell you why this colormap is grayscale? 

%% changing the color axis limits using set

% add a colorbar on the side
colorbar

colormap jet(64) % back to a normal colormap
% symmetric colorscales are the best
set(gca,'clim',[-30 30])
set(gca,'clim',[-3 3])
set(gca,'clim',[-1 1])



% mean-center the Gaussian so it has positive and negative values
gaus2d = gaus2d-mean(gaus2d(:));
imagesc(gaus2d)
colorbar

% asymmetric colorscales
set(gca,'clim',[0 1])
set(gca,'clim',[-1 0])

% watch what happens to the colorbar
colormap jet(64)
colormap jet(10)

%% saving figures as pictures

print(gcf,'-dpng','lots_of_red')


%% Do you want more introductory Matlab plotting code?
% Go to mikexcohen.com, download the material for the book
% "Analyzing Neural Time Series Data," and go through the script chapter04c.m
