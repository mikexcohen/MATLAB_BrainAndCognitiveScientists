function stats = mikeIsANiceGuy(data)
% mikeIsANiceGuy   Compute descriptive statistics on input data
%   stats = mikeIsANiceGuy(data) displays and returns the following descriptive statistics:
%         mean, median, minimum, maximum, count
%     INPUTS:
%          data   :   a 1D vector of numbers in single or double precision
%   
%     OUTPUTS:
%         stats   : a vector containing mean, median, minimum, maximum, count
% 
% code written by <your name here> <your email address here>

%% input checks go here...



%% compute the descriptive statistics

meanval   = mean(data);
medianval = median(data);
minval    = min(data);
maxval    = max(data);
N         = length(data);

% group outputs into one vector
stats = [ meanval medianval minval maxval N ];

%% display output stats for the user in order of output

disp([    'The mean value is ' num2str(meanval) '.' ])
disp([    'The median value is ' num2str(medianval) '.' ])
fprintf([ 'The minimum value is ' num2str(minval) '.' ]);
fprintf(  'The maximum value is %g.\n',maxval);
fprintf([ 'The number of numbers is ' num2str(N) '.\n' ])
disp(     'Done.')

%% end.
