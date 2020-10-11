%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 31
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% input patterns for backprop learning

% XOR patterns (first two numbers; the third number is for the bias term)
inputs = [ 1 0 1; 0 1 1; 1 1 1; 0 0 1 ];
output = [     0;     0;     1;     1 ];

% initialization
[outputs,predError] = deal( zeros(1,size(inputs,1)) );
totalError = nan(1,100);
figure(1), clf

%% setup model

nInputsNodes = size(inputs,2);
nHiddenNodes = 3;
nOutputNodes = 1;

% random initial weights
weights_i2h = randn(nInputsNodes,nHiddenNodes);
weights_h2o = randn(nHiddenNodes,nOutputNodes);

l_rate = .3;

%% run

toggle = true;
iteration = 0;
max_iterations = 300;

while toggle
    
    % loop through inputs
    for ini=1:size(inputs,1)
    
        
        %% forward part: compute inputs and errors
        
        % input-layer to hidden-layer: raw input -> weights -> sigmoid
        hdLayerResp = inputs(ini,:) * weights_i2h;
        hdLayerResp = 2./(1+exp(-hdLayerResp'*2))-1;
        
        
        % input-layer to output-layer: input -> weights -> sigmoid
        otLayerResp = hdLayerResp' * weights_h2o;
        otLayerResp = 2./(1+exp(-otLayerResp'*2))-1;
        
        
        % prediction error
        predError(ini) = otLayerResp - output(ini);
        
        
        % collect output-layer responses for plotting
        outputs(ini) = otLayerResp;
        
        %% backwards part: adjust weights based on error
        
        % adjust hidden -> output weights
        delta = l_rate * predError(ini) * hdLayerResp;
        weights_h2o = weights_h2o - delta;
        
        % adjust input -> hidden weights
        backprop = weights_h2o .* (1-hdLayerResp.^2) * inputs(ini,:);
        delta = l_rate * predError(ini) * backprop;
        weights_i2h = weights_i2h - delta';
        
    end
    
    iteration = iteration+1; % not auto-updated in a while loop!
    totalError(iteration) = sum(predError.^2);
    
    if totalError(iteration)<.01 || iteration>max_iterations
        toggle=false;
    end
    
    %% plot
    
    if mod(iteration,5)==0 % plot every 5th step
        subplot(221), cla
        plot(output,'bo','linew',2,'markersize',20,'markerfacecolor','b')
        hold on
        plot(outputs,'r*','linew',2,'markersize',10)
        set(gca,'xlim',[0 5],'ylim',[-.1 1.1],'xtick',1:4,'xticklabel',{'o1';'o2';'o3';'o4'})
        
        subplot(222)
        plot(totalError,'k','linew',2)
        set(gca,'xlim',[0 max_iterations],'ylim',[0 2.4])
        xlabel('Trials'), ylabel('Error^2')
        
        subplot(223)
        imagesc(weights_i2h), set(gca,'clim',[-2 2])
        title('input -> hidden weights')
        
        subplot(224)
        imagesc(weights_h2o'), set(gca,'clim',[-2 2])
        title('hidden -> output weights')
        
        drawnow
    end % end plotting
    
end % end simulation while-loop

%%


%% k-means clustering

load kmeans_data

% plot the raw data
figure(2), clf
plot(d(:,1),d(:,2),'ko')
axis([0 1 0 1])


% k-means clustering
k = 3; % how many clusters?
[groupidx,cents,sumdist,distances] = kmeans(d,k);
% for octave: kmeans is in the statistics package


figure(3), clf

% draw lines from each data point to the centroids of each cluster
lineColors = 'rkb';
hold on
for i=1:length(d)
    plot([ d(i,1) cents(groupidx(i),1) ],[ d(i,2) cents(groupidx(i),2) ],lineColors(groupidx(i)))
end

% now draw the raw data in different colors
for i=1:k
    plot(d(groupidx==i,1),d(groupidx==i,2),[ lineColors(i) 'o' ],'markerface','w')
end

% and now plot the centroid locations
% (why is this code after the previous plotting code?)
plot(cents(:,1),cents(:,2),'ko','markerface','g','markersize',10)


set(gca,'xlim',[0 1],'ylim',[0 1])
legend({'group 1';'group 2';'group 3'})

%% SVM in-sampling (!!)

% this and the next cell will not work in Octave. try the libsvm toolbox

load EEG_LR

% initialize
accu = zeros(size(timevec));
trueLabels = [ones(size(l_eeg,3),1); 2*ones(size(r_eeg,3),1)];


% loop over time (skip some points for speed)
for ti=1:20:size(l_eeg,2) 
    
    % organize the data for this time point
    data = squeeze( cat(3, l_eeg(:,ti,:), r_eeg(:,ti,:) ))';
    
    % fit the model
    svmModel = fitcsvm(data,trueLabels);
    
    % evaluate the model (note: training data used in testing, 
    % so this is in-sample testing)
    catLabel = predict(svmModel,data);
    
    % average accuracy
    accu(ti) = mean(catLabel==trueLabels);
end

% show a time course of the average accuracy results
figure(4), clf
plot(timevec(1:20:end),accu(1:20:end),'k','linew',2)
set(gca,'ylim',[.3 1])
hold on
plot(get(gca,'xlim'),[.5 .5],'k--')
xlabel('Time (ms)'), ylabel('Accuracy')

%% out-sample test

trueLabels = [ones(size(l_eeg,3),1); 2*ones(size(r_eeg,3),1)];
accu = zeros(length(timevec),length(trueLabels));

% warning! this takes a while to run...
for ti=1:20:size(l_eeg,2)
    
    % get data for this time point (doesn't change over trials!)
    data = squeeze( cat(3, l_eeg(:,ti,:), r_eeg(:,ti,:) ))';
    
    for triali=1:length(trueLabels)
    
        templabels = trueLabels;
        
        % remove test trial from training
        traindata = data;
        traindata(triali,:) = [];
        templabels(triali)  = [];
        
        % fit model on training data
        svmModel = fitcsvm(traindata,templabels);
        
        % Now predict only the trial that was not included
        % in the training set.
        catLabel = predict(svmModel,data(triali,:));
        accu(ti,triali) = catLabel==trueLabels(triali);
    end
end

% and plot
figure(4), clf
plot(timevec(1:20:end),mean(accu(1:20:end,:),2),'k','linew',6)

%% end.
