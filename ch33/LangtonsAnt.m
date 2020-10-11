
%% initialize

N = 500; % number of grid points

greed  = rand(N)<0;
curloc = [ round(N/20) round(N/3) ]; % starting location
antDir = 1; % numbers are somehow converted into directions...
nIters = 1000000; % number of iterations

vidspeed = 500;

%% setup figure

figure(1), clf
h = imagesc(rand(N));
axis xy, axis off, axis square, colormap gray

%% direction update matrices

dirupdates = [ 2  3  4  1 ;
               4  1  2  3 ];

locupdates = [ 1  0 -1  0 ;
               0 -1  0  1 ]';

%% run simulation

for ti=1:nIters
    
    % update ant's direction
    antDir = dirupdates( greed(curloc(1),curloc(2))+1 , antDir );
    
    % update current location
    curloc = curloc + locupdates(antDir,:);
    
    % wrap around edges
    curloc(curloc<1) = N;
    curloc(curloc>N) = 1;
    
    % update grid
    greed(curloc(1),curloc(2)) = 1-greed(curloc(1),curloc(2));
    
    % update screen
    if mod(ti,vidspeed)==1
        set(h,'CData',greed), drawnow
    end
    
end

%% end.
