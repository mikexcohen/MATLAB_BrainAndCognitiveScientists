%% MATLAB for Brain and Cognitive Scientists (MIT Press 2017)
% mikexcohen@gmail.com
%
% This script contains solutions and tips for most exercises in the book.
% If you are struggling with an exercise not in this list, email me with an
% explanation of what you've already tried and what you need help with.

%% Chapter 6

% exercise 10
((1:6)-1)*2+10

%% Chapter 7

% exercise 1
mod((1:15)-1,3)+1


% exercise 2
% Does each call to randn in the if-statement produce the same number?


% exercise 3
for i=1:20
    % define post-script
    ps = 'th'; % default
    if i==1
        ps = 'st';
    elseif i==2
        ps = 'nd';
    elseif i==3
        ps = 'rd';
    end
    
    % print text
    disp([ 'This is the ' num2str(i) ps ' element.' ])
end
% However, the above code is not a complete solution:
% If the loop goes above 100, it would print "101th element," which is
% incorrect. Thus, a full solution requires evaluating only the
% single-digit position of the number.


% exercise 5
for i=1:10
    % initialize
    if i==1
        data = zeros(10,40);
    end
    
    % rest of code here
end

%% Chapter 8

% exercise 2
a = randn(10,2);
a(4:end,:) = []

% exercise 6, hint:
[ 'the number ' num2str(3) '.' ]

%% Chapter 9

% exercise 3
figure(1), clf
h = patch(rand(4,1),rand(4,1),'k');
set(h,'FaceColor','y','EdgeColor','r','linewidth',5)
get(h) % for more options to change

%% Chapter 10

% exercise 6
A = [1 1; 2 2]; % for example...
Ainv = pinv(A);
A*Ainv % should be eye(2) for invertible matrix


% exercise 13
bsxfun(@times,randn(4,6),1:4)
% note: bsxfun is picky about vector orientation. 
%       In the code above, you need to transpose 1:4.

%% Chapter 11

% exercise 1
srate = 1000; % assume Hz
time  = 0:1/srate:10; % technically, this is *slightly* more than 10 seconds!
frex  = [ 2 8 11 29 ];
ampl  = [ 4 7  1 10 ];

% initialize
signal = zeros(size(time));
for fi=1:length(frex)
    signal = signal + ampl(fi)*sin(2*pi*frex(fi)*time);
end

% add noise (notice that the noise amplitude is a function of the average
% of all individual signal amplitudes)
signalLittle = signal + randn(size(time))*mean(ampl)*5;
signalLots   = signal + randn(size(time))*mean(ampl)*50;


figure(1), clf
subplot(311), plot(time,signal), title('Pure signal')
subplot(312), plot(time,signalLittle), title('Signal and little noise')
subplot(313), plot(time,signalLots), title('Signal and lots of noise')


% exercise 9
% Check out the coverart.m file!

% exercise 14
figure(1), clf
plot([0 2],[0 1],'k','linew',4)
axis([-3 3 -3 3])
axis square, grid on
xlabel('Real axis'), ylabel('Imaginary axis')

%% Chapter 12

% exercise 3

% import picture and convert to numeric
pic = imread('saturn.png');
pic = double(squeeze(pic(:,:,1)));

% create gaussian
[x,y]  = meshgrid(-250:250);
s = 30; % width, feel free to change
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

% making the images is up to you!

%% Chapter 13

% exercise 1

% 'measured' data
data      = round(randn(1,6)*5);
datatimes = 1:6;

% requested (interpolated) time points
newtimes  = 1:.001:6;

F = griddedInterpolant(datatimes,data,'spline');
newdata = F(newtimes);


%% Chapter 14

% exercise 1

% convert counts to probabilities
x  = .01:.01:.99;
dp = bsxfun(@minus,norminv(x)',norminv(x));
rb = -bsxfun(@plus,norminv(x)',norminv(x))/2;



% show the 2D d' space
figure(4), clf
imagesc(x,x,dp)
xlabel('False alarm rate')
ylabel('Hit rate')
hold on % we'll plot lines on top
axis square

% colors for isosensitivity curves
colorz = 'rbmk';

% the d' values
dp2plot = [ 1 1.5 2 2.5];
rb2plot = [.3 .5 .9 1.5];
tol = .01;

for dpi=1:length(dp2plot)
    
    % find points
    idx = find(dp>dp2plot(dpi)-tol & dp<dp2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(dp),idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',1,'markersize',6,'markerfacecolor',colorz(dpi))
    
    
    
    % find points
    idx = find(rb>rb2plot(dpi)-tol & rb<rb2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(rb),idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',1,'markersize',6,'markerfacecolor',colorz(dpi))
    
end

%% Chapter 15

% exercise 8

% generate the data
d = round(rand(10000,1)*(2345-14)+14);
% using the toolbox function
dd = prctile(d,98);

%% Chapter 16

% exercise 1

% inside the triple-loop will be the following
% create data (ns=number of points, cs=correlation strengths)
x = randn(ns,2);
x(:,2) = x(:,1)*cs(ci) + x(:,2)*sqrt(1-cs(ci)^2);

% correlation (or you can use corr or corrcoef
x = bsxfun(@minus,x,mean(x,1));
x = bsxfun(@rdivide,x,std(x,[],1));
tmp = x'*x / (ns(ni)-1);
corrmat(ni,ci,ti,1) = tmp(1,2); % store the absolute correlation (you'll also need the different from cs)


% exercise 6
a = 1:5; % row vector
a(:)


% exercise 7
realcor = corr(x,y);
for permi=1:1000
    permcor(permi) = corr(x,y(randperm(length(y))));
end
zcor = (realcor-mean(permcor)) / std(permcor);


%% Chapter 17

% exercise 6
% to compute eigenvector magnitude:
vecnorm = norm(evecs(:,10)) % for vector #10
% or 
vecnorm = sqrt(sum(evecs(:,10).^2))


%% Chapter 18

% exercise 6:
% The most common mistake here is to average together the complex Fourier
% coefficients over trials. You should first extract power with abs().^2 and
% then average the power values over trials.

% exercise 7:
data = randn(120,100);
dataHann = bsxfun(@times,data,.5*(1-cos(2*pi*(1:100)/(100-1))));
figure(5), clf
plot(sum( (dataHann-data).^2 ,1),'linew',2)
xlabel('Time (a.u.)'), ylabel('SSE')

% exercises 10:
load EEGrestingState.mat

n = length(eegdata);
epochs = reshape(eegdata,2048,[]);

fNot = mean( abs(fft(epochs)/2048).^2 ,2);
fHan = mean( abs(fft(bsxfun(@times,epochs,.5*(1-cos(2*pi*(1:2048)/(2048-1)))') )/2048).^2 ,2);

hz = linspace(0,srate/2,floor(2048/2)+1);

clf
plot(hz,fNot(1:length(hz)),hz,fHan(1:length(hz)),'linew',2)
set(gca,'xlim',[0 60])
xlabel('Frequency (Hz)'), ylabel('Power (\muV^2)')
legend({'Hann-tapered';'No taper'})


% exercise 11 (uses the same data as above)

fepoch = mean( abs(fft(epochs,n)/2048).^2 ,2);
fConti = abs(fft(eegdata)/n).^2;
hz = linspace(0,srate/2,floor(length(eegdata)/2)+1);


clf
plot(hz,fepoch(1:length(hz)),'linew',2), hold on
plot(hz,fConti(1:length(hz)),'linew',2)

fNot = mean( abs(fft(epochs)/2048).^2 ,2);
hz = linspace(0,srate/2,floor(2048/2)+1);

plot(hz,fNot(1:length(hz)),'linew',2)
legend({'epoched large N','continuous','epoched small N'})
set(gca,'xlim',[0 60])


%% Chapter 19

% Exercise 4

srate = 1000;
sigtime = 0:1/srate:2;
signal1 = sin(2*pi*30*sigtime);

wavtime = -2:1/srate:2;

% baseline time window
basetime = [.8 1.2];

nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

nFrex = 50;
frex  = linspace(1,srate/10,nFrex);
s     = linspace(4,12,nFrex) ./ (2*pi.*frex);
sigX  = fft(signal1,nConv);
tf    = zeros(nFrex,length(signal1));

for fi=1:nFrex
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    
    tp = abs( as(nHfkn:end-nHfkn+1) )*2;
    tf(fi,:) = 10*log10( tp/mean(tp(dsearchn(sigtime',basetime(1)):dsearchn(sigtime',basetime(2)))) );
end

figure(8), clf

subplot(211)
plot(sigtime,signal1)
set(gca,'ylim',[-1.05 1.05])

subplot(212)
contourf(sigtime,frex,tf,40,'linecolor','none')


% Exercise 9
% triangle frequency modulation
srate = 1000;
sigtime = 0:1/srate:8;
freqTS = abs(mod(sigtime,2)-1)*10;
meanF = mean(freqTS);
k = 2*pi/srate;
signal1 = sin(2*pi.*meanF.*sigtime + k*cumsum(freqTS-meanF));
% then do standard time-frequency analysis on signal1


















