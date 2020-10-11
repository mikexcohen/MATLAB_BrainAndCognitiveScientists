%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 6
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% examples of some basic uses of functions

% functions have outputs (x) and inputs (15,1)
% this function generates normally distributed random numbers.
x1 = randn(15,1);

% how does the previous output differ from this output?
x2 = randn(15);

% not limited to two dimensions!
x = randn(2,3,4)
whos x

% the function 'size' returns the number of elements along each dimension
% of the input (a.k.a., the size of the input)
sizeX = size(x);


% let's go back to a vector...
x = randn(15,1);

% some useful functions
minX = min(x);    % minimum
maxX = max(x);    % maximum
medX = median(x); % median value

%% using outputs as inputs

% is there a difference between this set of lines...
sizeX = size(x);
zmat  = zeros(sizeX);

% ... and this next line?
zmat = zeros( size(x) );
% above, I used the output of 'size' as an input into 'zeros'.


% Why doesn't this next line work? Does the error message help you understand?
zmat = zeros( max(x) );

% If you know why that line crashed, you should understand why this line
% doesn't (but why is this bad programming?).
zmat = zeros( round( max(x) ) );

% same line as previous but without the spaces. Which is easier to read?
zmat=zeros(round(max(x)));

%% functions with multiple inputs and multiple outputs

% multiple inputs (the second input is optional; is the first?)
size(x)
size(x,2)


% multiple outputs (the second one is optional; is the first?)
[val,idx] = max(x);


% max is an important function; make sure you understand its outputs!
[maxval,maxidx] = max([ 1 2 1 1 5 5 1 5 ])

% what output would you expect if you used 'min' instead of 'max'?

%% getting help in Matlab

% basic inline help
help mean

% usually more detailed help in its own window
doc mean % probably won't work in Octave


%% end.

