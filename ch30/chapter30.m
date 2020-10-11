%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 30
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% define timing parameters

srate = 10000; % sampling rate in Hz

% duration of simulation
sim_dur = 1; % in seconds
timevec = 0:1/srate:sim_dur - 1/srate;

%% define input stimulus

input = zeros(1,length(timevec));
input(dsearchn(timevec',.3):dsearchn(timevec',.7)) = 3;

%% define neuron properties

% voltage parameters
volt_rest   = -70; % resting membrane potential (mV)
volt_thresh = -50; % action potential threshold (mV)
volt_reset  = -75; % post-spike reset voltage (a bit below resting)

% membrane parameters
R_m = 10; % neuron membrane resistance (MOhm)
tau = 10; % time constant of decay (ms)

% initialize neuron membrane potential
neuronV = volt_rest + zeros(size(timevec));
spiketimes = [];

%% run the simulation

for timei = 1:length(timevec)-1
    
    %% test whether spike threshold was passed
    
    if neuronV(timei) > volt_thresh % threshold exceeded!
        % reset voltage
        neuronV(timei) = volt_reset;
        spiketimes = cat(1,spiketimes,timei);
    end
    
    %% update voltage in two steps
    
    % step 1: resting potential plus input scaled by membrane resistance
    restInput = volt_rest + input(timei)*R_m;
    
    % step 2: update voltage with peak
    neuronV(timei+1) = restInput + (neuronV(timei) - restInput) * exp(-1000/srate/tau);
    
end

% plotting niceties
neuronV(neuronV==volt_reset) = 40;

figure(1), clf
subplot(211)
plot(timevec,neuronV,'k','linew',2)
set(gca,'ylim',[-100 50])
xlabel('Time (s)'), ylabel('Voltage (mV)')

subplot(212)
plot(timevec,input)
set(gca,'ylim',[min(input)-.1 max(input)*1.1])
xlabel('Time (s)'), ylabel('Input current')

%% and now for a population...

%% (re)define timing parameters

srate = 10000; % sampling rate in Hz

N_exc = 80; % number of excitatory neurons
N_inh = 20; % number of inhibitory neurons

% duration of simulation
sim_dur = 1; % in seconds
timevec = 0:1/srate:sim_dur - 1/srate;

%% define input stimulus

input = zeros(1,length(timevec));
input(dsearchn(timevec',.3):dsearchn(timevec',.7)) = 3;

%% define neuron properties

% voltage parameters
volt_rest   = -70; % resting membrane potential (mV)
volt_thresh = -50; % action potential threshold (mV)
volt_reset  = -75; % post-spike reset voltage (a bit below resting)

% membrane parameters
R_m = 10; % neuron membrane resistance (MOhm)
tau = 10; % time constant of decay (ms)

% initialize neuron membrane potential
neuronV = volt_rest + zeros(N_exc+N_inh,length(timevec));

%% run the simulation

for timei = 1:length(timevec)-1
    
    %% test whether spike threshold was passed
    
    spikedNeurons = neuronV(:,timei) > volt_thresh; % threshold exceeded!
    
    % reset voltage
    neuronV(spikedNeurons,timei) = volt_reset;
    %spiketimes = cat(1,spiketimes,timei);
    
    %% update voltage in two steps
    
    % step 1: resting potential plus input scaled by membrane resistance
    restInput = volt_rest + input(timei)*R_m;
    restInput = [restInput*ones(1,(round(N_exc+N_inh)/2)) volt_rest*ones(1,(round(N_exc+N_inh)/2))]';
    
    % step 1.5: integrate output of other neurons
    restInput = restInput + 25*sum(spikedNeurons(1:N_exc)) - 10*sum(spikedNeurons(N_exc+1:end));
    
    % step 2: update voltage with peak
    neuronV(:,timei+1) = restInput + (neuronV(:,timei) - restInput) * exp(-1000/srate/tau);
    
end

%% and plot

figure(2), clf

% an image of the membrane potentials of all neurons
subplot(211)
imagesc(timevec,[],-(neuronV-mean(neuronV(:))))
xlabel('Time (s)'), ylabel('Neuron number')
colorbar

% plotting niceties
neuronV(neuronV==volt_reset) = 40;

subplot(212)
h = plot(timevec,neuronV([4 94],:));
xlabel('Time (s)'), ylabel('Voltage (mV)')
set(gca,'ylim',[-80 50])

% question: why do the colors in the top plot change 
% when you re-run this cell?

%% Izhikevich neurons 
% code adapted from Izhikevich 2003

% parameters that control neuron's behavior
a = .03;
b = .25;
c = -60;
d =   4;

tau = .25; % what is this value in Hz?
tspan = 0:tau:1000;

% define time series of input
T1 = zeros(size(tspan));
T1(dsearchn(tspan',200):dsearchn(tspan',800)) = 1;

V = -70;
u = b*V;
[VV,uu] = deal(zeros(size(T1)));

for ti=1:length(tspan)
    
    % membrane potential
    V = V + tau*(.04*V^2 + 5*V + 140 - u + T1(ti) ); % T1 is the input
    u = u + tau*a*(b*V-u);
    
    if V > 30 % there was a spike
        VV(ti+1)=30;
        V = c;
        u = u + d;
    else % there was no spike
        VV(ti+1)=V;
    end
    uu(ti+1)=u;
end

% and plot it
figure(3), clf
plot(tspan,VV(1:end-1),tspan,10*T1-88);
axis([0 max(tspan) -90 35])
xlabel('Time (ms)'), ylabel('Membrane potential')

%% Rescorla-Wagner-esque learning model

nTrials = 100;
lrate = .3;

rewProbs = [.7 .2];

w = .5+zeros(nTrials+1,2);
[action,rewpred,pPickAct1] = deal( zeros(1,nTrials) );

for triali=1:nTrials
    
    % compute probability of picking action
    pPickAct1(triali) = exp(w(triali,1)) / sum(exp(w(triali,:)));
    
    % pick an action based on weighted probability
    action(triali) = 1 + (pPickAct1(triali)<rand);
    
    % is this action rewarded? (convert to -1/+1)
    reward = rand < rewProbs(action(triali));
    
    % compute prediction error (aka delta)
    rewpred(triali) = reward - w(triali,action(triali));
    
    % update weights for the next trial
    w(triali+1,action(triali)) = w(triali,action(triali)) + lrate*rewpred(triali);
    w(triali+1,3-action(triali)) = w(triali,3-action(triali));
end

figure(4), clf
subplot(311), plot(w,'linew',2), legend({'w1';'w2'})
xlabel('Trial'), ylabel('Action weights')
set(gca,'xlim',[1 nTrials])

subplot(312), plot(1:nTrials,rewpred,'linew',2)
hold on, plot(get(gca,'xlim'),[0 0],'k:','linew',3)
xlabel('Trial'), ylabel('Reward prediction error')
set(gca,'xlim',[1 nTrials])

subplot(313), plot(1:nTrials,pPickAct1,'linew',2), set(gca,'ylim',[0 1])
hold on, plot(get(gca,'xlim'),[.5 .5],'k:','linew',3)
xlabel('Trial'), ylabel('Action probability')
set(gca,'xlim',[1 nTrials])

%% end
