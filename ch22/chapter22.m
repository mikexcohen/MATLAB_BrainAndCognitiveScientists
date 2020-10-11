%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 22
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% analog vs. binary representation of a spike

% use data from chapter 23
load ../ch23/times_090425blk10_ch115.mat

figure(1), clf

% analog
subplot(211)
plot(spikes(10,:),'k','linew',2)
set(gca,'xlim',[0 13])
xlabel('Time (a.u.)')
ylabel('Voltage (mV)')


% binary simplification
subplot(212)
set(gca,'xlim',[0 13],'ylim',[-.1 1.2],'ytick',[0 1],'yticklabel',{'no';'yes'})
hold on
plot([1 12],[0 0],'k','linew',2)
plot([5 5],[0 1],'k','linew',2)
xlabel('Time (a.u.)')
ylabel('Spike?')

%% import data

% these data were downloaded from:
% https://crcns.org/data-sets/motor-cortex/alm-1/about-alm-1

fid = fopen('03-04-11-aa_sig1_spikes.dat');

spikenum = 1;
trialnum = 0;
APsS = zeros(100,2);
APsM = zeros(1,1001);

while ~feof(fid)
    
    % get a line of code
    d=fgetl(fid);
    
    % if it's an empty line, skip to the next line
    if isempty(d), continue; end
    
    % if this is an inter-trial comment, increment trial number 
    % and then skip forward
    if d(1)=='#'
        if strcmpi(d(1:5),'# int')
            trialnum = trialnum+1;
        end
        continue;
        
        
    % otherwise, it's real data; grab the spike time 
    % and update the matrices
    else
        spiketime = ceil(sscanf(d,'%g')/1000);
        
        APsS(spikenum,:) = [trialnum spiketime];
        APsM(trialnum,spiketime) = 1;
        spikenum = spikenum+1;
    end
end
fclose(fid);

npnts = size(APsM,2);
timevec = 1:npnts;


% plot the data as an image
figure(2), clf
subplot(311)
imagesc(bsxfun(@times,APsS',[10 1]')), colormap gray
xlabel('Spike count')
set(gca,'ytick',1:2,'yticklabel',{'Trial';'Time'})

% and a different kind of image
subplot(312)
plot(APsS(:,2),APsS(:,1),'k.','markersize',3)
set(gca,'xlim',timevec([1 end]),'ylim',[1 trialnum])
xlabel('Time (ms)'), ylabel('Trial')

% ...the same kind of image but using a different format
subplot(313)
imagesc(timevec,[],2-APsM)
xlabel('Time (ms)'), ylabel('Trial')

%% spike rate (spikes/second)

tidx = dsearchn(timevec',[0 1000]');

dt = timevec(tidx(2)-tidx(1)) / 1000;
spikeRates = sum(APsM(:,tidx(1):tidx(2)),2) / dt;

figure(3), clf
plot(spikeRates)
xlabel('Trial'), ylabel('spike rate (sp/s)')
set(gca,'xlim',[0 size(APsM,1)+1])

%% PSTH

figure(4), clf

% PSTH at the temporal resolution of the data
subplot(211)
dt = mean(diff(timevec)) / 1000;
plot(timevec,mean(APsM)./dt,'k','linew',2)
set(gca,'xlim',timevec([1 end]))
xlabel('Time (ms)'), ylabel('Spike rate')


% downsampling to a smaller number of bins
subplot(212)
bins = ceil(timevec/5);
[spPerbin,timevecBin] = deal( zeros(1,max(bins)) );
for i=1:length(spPerbin)
    spPerbin(i) = mean( sum(APsM(:,bins==i),2) ,1);
    timevecBin(i) = mean(timevec(bins==i));
end

dt = mean(diff(timevecBin))/1000;
plot(timevecBin,spPerbin./dt,'k','linew',2)
set(gca,'xlim',timevecBin([1 end]))
xlabel('Time (ms)'), ylabel('Spike rate')

%% end.
