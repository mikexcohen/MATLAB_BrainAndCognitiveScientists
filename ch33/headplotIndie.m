% headplot() - plot a spherically-splined EEG field map on a semi-realistic
%              3-D head model. Can 3-D rotate the head image using the left
%              mouse button.

%
% Authors: Arnaud Delorme, Colin Humphries, Scott Makeig, SCCN/INC/UCSD,
%          La Jolla, 1998-

% this function was modified by mikexcohen@gmail.com.
% Minor modifications were made to allow this file to be independent of the eeglab
% toolbox. All credit for the function belongs with eeglab!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Arnaud Delorme, Colin Humphries and Scott Makeig,
%               CNL / Salk Institute, Feb. 1998
%
% Spherical spline method: Perrin et al. (1989) Electroenceph clin Neurophys
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function headplotIndie(values, spline_file, clim)

%%%%%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULT_LIGHTS = [-125  125  80; ...
    125  125  80; ...
    125 -125 125; ...
    -125 -125 125];    % default lights at four corners

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open head mesh and electrode spline files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(spline_file, '-mat');

values = values(indices);

% load mesh file
% --------------
load('mheadnew.mat');
index1 = sort(unique(TRI1(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%

meanval = mean(values); values = values - meanval; % make mean zero
lamd = 0.1;
C = pinv([(G + lamd);ones(1,length(values))]) * [values(:);0]; % fixing division error
P = zeros(1,size(gx,1));
for j = 1:size(gx,1)
    P(j) = dot(C,gx(j,:));
end
P = P + meanval;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%
cla % clear axis

W = zeros(1,size(POS,1));
m = 64;


idx = min(m,round((m-1)*(P-clim(1))/(clim(2)-clim(1)))+1); % get colormap indices

W(index1) = idx;
colormap(jet(m))
p1 = patch('Vertices',POS,'Faces',TRI1,'FaceVertexCdata',W(:),...
    'FaceColor','interp', 'cdatamapping', 'direct', 'tag', 'mesh');    %%%%%%%%% Plot scalp map %%%%%%%%%

axis([-125 125 -125 125 -125 125])
axis off % hide axis frame

%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn on lights
%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(DEFAULT_LIGHTS,1)
    light('Position',DEFAULT_LIGHTS(i,:),'Color',[1 1 1],'Style','infinite');
end

set(p1,'DiffuseStrength',.6,'SpecularStrength',0,'AmbientStrength',.3,'SpecularExponent',5,'EdgeColor','none','vertexnormals', NORM)
lighting phong  % all this gives a matte reflectance

rotate3d on
axis image    % keep the head proportions human and as large as possible

