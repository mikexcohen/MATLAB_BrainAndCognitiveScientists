%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 7
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% if-then etc

% basic structure of if-then statements
if 2>1
    disp('I like big numbers.')
end

%% illustrating elseif and else components of if-then statements

% see book for explanation of letters

r = randn; % with no inputs, a single random number
if r>0 % A
    disp([ num2str(r) ' is positive.' ])
elseif r>2 % B
    disp('Snuffleupagus is real!')
elseif r<0 && r>-1 % C
    disp('r is small but I''ve seen smaller');
else % D
    disp('r is really small')
end

%% nested if statements

% if statements can appear inside if statements

r = rand;
if r>0
    disp('r is definitely positive')
    if r>2, disp('>2'); end
    if r>1
        disp('>1'); end % don't do this!
else
   disp('r is most likely negative')
end

%% for-loops

% count from 1 to 10 in integer steps...
for i=1:10
    % and print out the number on each iteration.
    disp(num2str(i))
end

% another example, using non-integer counting
for i=-5:1.1:6
   disp([ 'My new favorite number is ' num2str(i) '.' ])
end


% the looping variable can be used as an index and as a variable
for i=1:5
    sq(i,1) = i;
    sq(i,2) = i^2;
end
sq

%% tip for naming looping variables

% good practice is to have the variable name end in i (i=index)
% Notice the looping index variable names. And notice the informative ends.

% loop over channels
for chani=1:nChannels
    
    % <some content here>
    
    % loop over conditions
    for condi=1:nConditions
        
        % <dozens of lines of content here>
        
        % loop over frequencies
        for freqi=1:nFrequencies
            
            % <yet more content here>
            
            % enter data into matrix
            resmatrix(chani,condi,freqi) = analysis_results;
            
        end % end frequency loop
        
        % <content continues>
        
        % more data into another matrix
        moreMatrix(chani,condi) = otherstuff;
        
    end % end condition loop
end % end channel loop

%% skipping an iteration in a loop

for i=1:10
    if mod(i,2)==0, continue; end
    disp([ 'Running iteration ' num2str(i) '.' ])
end

% This illustrates how the 'mod' function works.
% The modulus is the remainder after division.
[ 1:10; mod(1:10,2) ]

%% while-loops

% a never-ending story:
i=1;
while true
    i=i+1;
    disp(num2str(i))
end % in Octave, nothing will print until after breaking the loop


% this loop is for people with sub-infinite patience
i=1;
while i<50
    i=i+1;
    disp(num2str(i))
end

%% breaking out of while-loops

% The difference between this and the following loop is subtle.
% Note the value of 'r' at the end of each loop.

% this one uses a toggle to stop the loop
r = 75;
toggle = true;
while toggle
    r = r/1.1;
    if r<1, toggle=false; end
    disp([ 'Current value is ' num2str(r) ])
end




% this one uses a hard break
r = 75;
while true
    r = r/1.1;
    if r<1, break; end
    disp([ 'Current value is ' num2str(r) ])
end

%% try-catch statements

% define a 2-element vector
e = [1 4];

% this produces an error (of course)
e(3)>6

% now watch:
try
    if e(3)>6, disp('e(3) is big.'), end
catch me;
end
% No error. What is me (the variable)?


% if-then statement nested inside a try-catch statement
try
    if e(3)>6, disp('e(3) is big.'), end
catch me;
    if e(2)>6, disp('e(2) is big.')
    else disp('e(1) is my friend.')
    end
end

%% end.
