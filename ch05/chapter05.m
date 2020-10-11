%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 5
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% starter introduction to variables

% a variable is a place-holder for information. 
% To create a variable, assign to it a value using the = sign:
asdf = 7

% Matlab repeats what you wrote. Try adding a semi-colon at the end:
asdf = 7;

% By the way, this is a comment. Comments are lines that start with a percent sign.
% Matlab ignores comments, so you can use them to make notes about what the code does.
asdf = 8; % comment after code

%% whos

% type this line in the command window:
whos

% then this line:
clear asdf

% and then whos again:
whos

%% variable names

% illegal!
4asdf = 4;
#asdf = 4; % this one's OK in Octave because # can also be used to indicate comments

% legal and useful:
data2lookAt = 8; % number in the middle of a variable
spike_results = 150; % underscore to separate words

% legal but potentially confusing:
asdf = 4; % lower-case s
aSdf = 5; % upper-case S

% legal but also not great choices:
a = 4;
i = 1;

% also legal, but a terrible idea:
mean = 4; % 'mean' is a Matlab function!

%% arrays and matrices

% you can clear variables during a script
clear a i 
clear a* % the 'wildcard' * clears all variables that start with 'a'

anArray = [1 3 2 5 3 6 7 9];

aMatrix = [1 3 2; 5 3 6; 7 8 9];

whos

%% basic mathematical expressions

% create variables
a = 4;
b = 3;

% Perform simple mathematical operations using those variables.
% Run each line at a time.
a+b
a*b
(a*4) + (b/2)
a*4+b/2

%% boolean (true/false) variables

isEverythingOK = 1==2

% this variable contains 'FALSE' although it looks like a zero
whos isEverythingOK % note the "Class" of the variable.

% Matlab can convert logicals to numbers on the fly. 
% But this can cause confusion; use this feature sparingly.
var = true;
var2 = var + 4

whos var* % again, note the variable classes

%% string variables

% Non-alphanumeric characters are fine here, because they are in
% the string (the content of the variable), not in the variable name.
nameVariable = 'Mike X Cohen@#$%';

% potential confusing resulting from Matlab storing strings 
string4 = '4';
number4 = 4;
number4 + string4

%% cells

% square brackets can be used to concatenate, but the result is probably
% not what you hoped it would be in this case:
a = [52 'hello' 52]

% instead, a cell array will keep each piece of information separate.
celery{1} = 52;
celery{2} = 'hello';
celery{3} = 52;

% the function 'cell' will initialize an empty cell matrix, in this case,
% 2 rows and 4 columns
celery2d = cell(2,4)

%% structures

% a 'structure' is a variable that can store different kinds of information
mouse.name          = 'mickey';
mouse.surgeryDate   = '24 March';
mouse.geneType      = 'wildtype';
mouse.age           = 65; % in days since birth
mouse.numElectrodes = 3;


% and now we add a second element to the structure
% Note the difference in spacing here and above;
% spacing isn't necessary but can improve visibility.
mouse(2).name = 'minney';
mouse(2).surgeryDate = '2 April';
mouse(2).geneType = 'GAD2Cre';
mouse(2).age = 82; % in days since birth
mouse(2).numElectrodes = 6;



% access all data from one element of the structure:
mouse(1)

% or access information from one field from all elements:
mouse.age
[mouse.age]
{mouse.age}
% notice how Matlab differently concatenates the previous three lines

%% the colon operator

% watch Matlab skip
1:5  % same thing has 1:1:5

% and skip by 2's
1:2:5


% use square brackets to skip some elements
[1:3 6 8:13 12:-1:5]

% we can use these numbers to create a variable
indices = [ -3:.1:-2.7 15:39:200 ];

%% using the colon operator for indexing

numbers = 1:2:200;

% access the 5th element
numbers(5)

% access elements 3-5
numbers(3:5)

% combine square-bracket-concatenation with the color operator
numbers([1:3 6 8:13])

% why does this line crash?
numbers(2.5:.02:3)


% and now combine logical variables with indexing (a.k.a. logical indexing)
logicalIdx = 1:100>50;

numbers(logicalIdx)


% uh oh... why doesn't this work?
numbers(double(logicalIdx))


% doesn't crash, but gives strange results. Why?
logicalIdx = 1:100>0;
numbers(double(logicalIdx))

numbers(double(logicalIdx)+1)


%% initializing matrices

% is there a difference between these next two lines?
mat1 = zeros(3);
mat2 = zeros(3,3);


% other ways to initialize
mat3 = ones(2,4);
mat4 = 7.43*ones(2,4) + 10;
mat5 = true(1,3);
mat6 = nan(8,2,3,4,5,2);

%% soft-coding

% The variable nbins is 'soft-coded' because it is defined once and then
% used throughout the rest of the script.
nbins = 10;

% <lots of code...>
results_mat = nan(nbins,1);

% <lots more code...>

% a loop... you'll learn about this in Chapter 7
for bini=1:nbins
   %<even more code...>
end

%% now the same script is hard-coded

% <lots of code...>
results_mat = nan(10,1);

% <lots more code...>

for bini=1:10
   %<even more code...>
end

%% Do you want more introductory Matlab code?
% Go to mikexcohen.com, download the material for the book 
% "Analyzing Neural Time Series Data," and go through the script chapter04a.m

