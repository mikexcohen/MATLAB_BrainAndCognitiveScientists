%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 28
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% average as a least-squares problem

n = 10; % number of data points

% generate some data
b = linspace(1,3,10)' + rand(n,1);

% the "design matrix" is just a vector of ones
A = ones(n,1);

% solution to the least squares problem
x = (A'*A)\A'*b;

% check that least-squares and mean are the same
fprintf('\n  Least squares solution: %g \n  Matlab ''mean'' function: %g \n\n',x, mean(b))

%% Improper design matrix for two independent variables

A = ones(n,2);
x = (A'*A)\A'*b;

%% The correct design matrix for an intercept and linear effect

A = [ ones(n,1) (1:n)' ];
x = (A'*A)\A'*b;

% why do these next two lines give errors?
A = [ ones(10,1) 1:n ];
A = [ ones(10,1) 1:n' ];

%% compute the model-predicted data 

% yHat are the data predicted by the model 
% (design matrix scaled by the coefficients)
yHat = x(1)*A(:,1) + x(2)*A(:,2);
% Note: in this case, A(:,1) is all ones, 
% so the multiplication here is not necessary.


figure(1), clf

% important visual check: plot the original data and model-predicted data
subplot(211)
plot(1:n,b,'o','markerface','b','markersize',20)
hold on
plot(1:n,yHat,'rp-','linew',2,'markersize',30,'markerface','k')
title('Correct')

set(gca,'xlim',[0 n+1],'ylim',[1 4.5])
legend({'observed data';'predicted data'})

%% two common mistakes to avoid

% Swapping the coefficients and design matrix columns.
yHat_error1 = x(2)*A(:,1) + x(1)*A(:,2);

% coefficients scale the y-axis data, not the x-axis data
yHat_error2 = x(1)*b + x(2)*b;

subplot(223)
plot(A(:,2),b,'o','markerface','b','markersize',20)
hold on
plot(A(:,2),yHat_error1,'rp-','linew',2,'markersize',30,'markerface','k')
set(gca,'xlim',[0 n+1],'ylim',[1 max(yHat_error1)+1])
title('Error #1')

subplot(224)
plot(A(:,2),b,'o','markerface','b','markersize',20)
hold on
plot(A(:,2),yHat_error2,'rp-','linew',2,'markersize',30,'markerface','k')
set(gca,'xlim',[0 n+1],'ylim',[1 max(yHat_error1)+1])
title('Error #2')

%% R-squared

resvar = sum( (b-yHat).^2 );
totvar = sum( (b-mean(b)).^2 );

r2 = 1 - ( resvar / totvar );

subplot(211)
title([ 'R^2 of model to data: ' num2str(r2) ])

%% Polynomial fitting

iv = (1:n)'; % iv = independent variable
regcoefs = polyfit(iv,b,1);

% check that least-squares and mean are the same
fprintf('\n    ''Manual'' least squares: [%g, %g] \n Matlab ''polyfit'' function: [%g, %g] \n\n',x, regcoefs)

%% More polynomial fitting

srate = 1000;
time = 0:1/srate:6;
F = griddedInterpolant(0:6,100*rand(7,1),'spline');
data = F(time) + 10*randn(size(time));

% polynomial fit
polycoefs = polyfit(time,data,7);

% estimated data based on the coefficients
yHat = polyval(polycoefs,time);


figure(2), clf
plot(time,data,'k'), hold on
plot(time,yHat,'r','linew',5)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')

legend({'Original data';'Predicted data'})

%% Using polynomial fitting to isolate slow vs. fast signal components

figure(3), clf

subplot(211)
plot(time,yHat,'k')
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Isolating the slow part of the signal')

subplot(212)
plot(time,data-yHat,'k')
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Isolating the fast part of the signal')

%% load EEG data and extract reaction times in ms

load /home/mxc/Dropbox/books/ants/code/sampleEEGdata.mat

rts = zeros(size(EEG.epoch));

% loop over trials
for ei=1:EEG.trials
    
    % find the index corresponding to time=0, i.e., trial onset
    [~,zeroloc] = min(abs( cell2mat(EEG.epoch(ei).eventlatency) ));
    
    % reaction time is the event after the trial onset
    rts(ei) = EEG.epoch(ei).eventlatency{zeroloc+1};
end

% create design matrix
A = [ ones(EEG.trials,1) rts' ];

%% define convolution parameters for time-frequency analysis

freqrange  = [2 20]; % extract only these frequencies (in Hz)
numfrex    = 30;     % number of frequencies between lowest and highest


% set up convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;
nCyc    = logspace(log10(4),log10(12),numfrex);

% create wavelets
cmwX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    s   = nCyc(fi) / (2*pi*frex(fi));
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / (2*s.^2) );
    
    % compute fourier coefficients of wavelet and normalize
    cmwX(fi,:) = fft(cmw,nConv);
    cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
end


% initialize time-frequency output matrix
tf = zeros(numfrex,EEG.pnts);
tf3d = zeros(numfrex,EEG.pnts,EEG.trials);

% compute Fourier coefficients of EEG data (doesn't change over frequency!)
eegX = fft( reshape(EEG.data(47,:,:),1,[]) ,nConv);

% loop over frequencies
for fi=1:numfrex
    
    % second and third steps of convolution
    as = ifft( cmwX(fi,:).*eegX ,nConv );
    
    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % extract power from all trials
    tf3d(fi,:,:) = abs(as).^2;
    
end % end frequency loop

%% now compute correlations

% reshape the 3D matrix to 2D
tf2d = reshape(tf3d,numfrex*EEG.pnts,EEG.trials)';

% the 2D matrix can be used in a single least squares equation
x = (A'*A)\A'*tf2d;
covmat = reshape(x(2,:),numfrex,EEG.pnts); % demushing

%% show the design and data matrices

figure(4), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(A)
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials')


ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d)
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials')
xlabel('Timefrequency')

colormap gray

%% show the results

figure(5), clf

% show time-frequency map of regressors
contourf(EEG.times,frex,covmat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-200 1200],'clim',[-.012 .012])

%% transforms to linearize nonlinear data

n = 10000;


% to use this cell, uncomment one of the following code sections

% %%% squared to  uniform via sqrt
% orig_data  = linspace(1,100,n).^2 + 10*rand(1,n);
% trans_data = sqrt( orig_data );
% 
% %%% uniform (0 to 1) to normal via atanh
% orig_data  = rand(1,n);
% trans_data = atanh( 2*orig_data-1 ); % atanh is a.k.a. Fisher-Z transform
% 
% %%% powers to normal via log
% orig_data  = 2.^randn(1,n); % try with other integers
% trans_data = log( orig_data );
% 
%%% anything to uniform (0 to 1) via rank
orig_data  = 2.^randn(1,n); % try with other integers
trans_data = tiedrank( orig_data ); % tiedrank in stats toolbox or (for octave) statistics package
% and then uniform to normal via atanh
trans_data = atanh( 2*(trans_data./max(trans_data))-1);



% and now for plotting5
figure(6), clf

subplot(131)
hist(orig_data,50)
set(gca,'xlim',[min(orig_data) max(orig_data)])
xlabel('Values'), ylabel('Count')
title('Original values')

subplot(132)
hist( trans_data ,50)
set(gca,'xlim',[min(trans_data) max(trans_data)])
xlabel('Values'), ylabel('Count')
title('Transformed values')

subplot(133)
plot(orig_data,trans_data,'ko')
xlabel('Original values'), ylabel('Transformed values')

%% end.
