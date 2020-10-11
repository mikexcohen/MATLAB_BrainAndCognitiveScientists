%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 8
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% loading in .mat files

% to load this file, Matlab needs to be in the same directory as the file,
% or you need to have the folder where the file lives in Matlab's path.
load ch08filename.mat

% full path specification also works:
%load /home/mxc/ch08filename.mat % for unix
%load C:/Users/mxc/Desktop/ch08filename.mat % for windows

% whenever you import a new file, always see what's contained in the file!
whos

% alternative method
clear
datafile = load('ch08filename.mat')

%% loading multiple files into a cell array

% initialize!
data = cell(15,1)

% 15 files, all with the same filename pattern
for filei=1:15
    data{filei} = load([ 'data_rat' num2str(filei) '.mat' ]);
end

%% getting information about folder contents

allFiles  = dir; % nonspecific
someFiles = dir('*rat*.mat'); % specific to this pattern

% you can also use paths in the dir command
files = dir('/path/ch08/*rat*.mat');

% ... although the output of dir does not contain the full path
files(3).name

% This is a great opportunity to apply soft-coding!
% list directory and get files
filedir = '/home/mxc/mbcs/ch08/';
files   = dir([ filedir '*rat*.mat' ]);

% initialize
data  = cell(size(files));

% loop over files and import one-by-one
for filei=1:length(files)
    data{filei} = load([ filedir files(filei).name ]);
end

%% simple graphical user interface (GUI) to select files

[filename,filepath] = uigetfile('*rat*.mat');

% similar function to select just a folder
filepath = uigetdir;

%% saving data in a .mat file

% Why does the next line crash? How can you make it work?
save('output_filename.mat','var1','var2','var3');

%% an example loop used to import data, run analyses, and save results.

% (Note: this code is for viewing purposes; it's not designed to work.)

% get file names
datadir = '/path/mbcs/ch08/';
filz = dir([ datadir '*_file2analyze.mat' ])

% loop over files
for filei=1:length(filz)

	% specify results output file name
    outfilename = [datadir filz(filei).name(1:end-16) 'results.mat' ];

    % if that file exists, the analysis has been run, so don't re-run
    if exist(outfilename,'file'), continue; end

    % load in the data
    load([ datadir filz(filei).name ])

    % perhaps hundreds of lines of code here...
    % <analyze...>

    % and save the results to a file
    save(outfilename,'var1','etc')
end % end loop over files

%% importing simple text files

% open this text file in a text editor and compare it to 'data'
data = load('textdata_easy.txt');

%% advanced importing text data

% Here we borrow from C language to flexibly read in mixed data. Let's say
% you have some poorly organized behavioral data files to read in, but at 
% least you know what text strings to look for: 

fid = fopen('headache_data.txt','r');
% fid is a pointer to a location on the physical hard disk (similar to how
% we used variables as handles to axes when plotting). The 'r' means read
% (later we'll use 'w' for write).

% In this particular example, we will extract the trial number, subject
% choice, reaction time (RT), and accuracy for each trial. Fields are separated by tabs.

behavioral_data=[]; % initialize... we can't initialize the full matrix, because we don't know how big this will be.

% The following code will remain inside a loop, reading in and processing new
% lines of data, until we reach the end of the file.
datarow=1;

while ~feof(fid) % feof tests whether we're at the end of the file.
    
    dataline = fgetl(fid); % read a line ("file get line")
    
    dataline = regexp(dataline,'\t','split');
    % regexp can be used to cut data according to delimiters. Here we will
    % cut this string of characters into a cell array in which elements of
    % the array are separated by tabs.
    % See next cell for more examples using regexp.
    
    % here we use strcmpi to compare strings. The "i" means to ignore case.
    if ~any(strcmpi('trial',dataline))
        continue % continue means to skip to the next iteration of the loop.
    end
    
    trial_column    = find(strcmpi('trial',   dataline));
    choice_column   = find(strcmpi('choice',  dataline));
    rt_column       = find(strcmpi('rt',      dataline));
    accuracy_column = find(strcmpi('accuracy',dataline));
    
    behavioral_data(datarow,1) = str2double(dataline{trial_column+1});      % Note that we didn't initialize the size of the variable "behavioral_data" so matlab gives a warning.
    behavioral_data(datarow,2) = str2double(dataline{choice_column+1});     % If the variable is relatively small, it doesn't matter. 
    behavioral_data(datarow,3) = str2double(dataline{rt_column+1});         % If the variable is large, however, it's best to initialize it to something really big, and then cut it down to size afterwards.
    behavioral_data(datarow,4) = str2double(dataline{accuracy_column+1});   % See chapter 4 in the book for further discussion of matrix initializations.
    
    datarow=datarow+1; % increment row
end

fclose(fid); % don't forget to close the file after you finish it!

%% regexp (regular expression)

sentence = 'Hello my name is Mike and I like 8 bikes.';
s1 = regexp(sentence,' ','split') % parse by spaces
s2 = regexp(sentence,'M','split') % parse by M

%% converting strings to numbers

% the goal here is to convert the string '42' into the number 42

aNumber = '42';
str2double(aNumber)
sscanf(aNumber,'%g')

% but these methods won't work on cell arrays:
str2double(aNumber)
sscanf(aNumber,'%g')

% instead, either use a loop or the function cellfun
% cellfun means a function (indicated by @) is applied to each cell
cellfun(@str2double,s1);

%% exporting text files: the easy way

dlmwrite('filename',data,'\t');
% What is the extension of the file on your computer?

%% exporting text files: the flexible way

% file ID pointer
fid = fopen('data2write.txt','w'); % 'w' instead of 'r'!

% print column headers in the first line
fprintf(fid,'accuracy\treactionTime\n');

% loop through rows in the matrix
for ri=1:size(data,1)
	% and write each row of data to the text file
    fprintf(fid,'v1: %g\tv2: %g\n', data(ri,1),data(ri,2) );
end

% close the file when finished!
fclose(fid);

%% importing images

% import a picture file and store as a matrix
% this image comes with Matlab, and is provided in the download folder for
% Octave users
pic = imread('saturn.png'); % comes with Matlab

% show the contents of the image
% you'll learn more about plotting in the next chapter
subplot(221)
% plot the whole image
imagesc(pic)

% note for Octave users: depending on your version of Octave, you might
% need to use imshow instead of imagesc

% loop through each RGB layer and show that slice of the image
for i=1:3
    subplot(2,2,i+1)
    imagesc(squeeze(pic(:,:,i)))
end

%% end.
