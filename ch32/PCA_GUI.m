function varargout = PCA_GUI(varargin)
% PCA_GUI MATLAB code for PCA_GUI.fig
%      PCA_GUI, by itself, creates a new PCA_GUI or raises the existing
%      singleton*.
%
%      H = PCA_GUI returns the handle to a new PCA_GUI or the handle to
%      the existing singleton*.
%
%      PCA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCA_GUI.M with the given input arguments.
%
%      PCA_GUI('Property','Value',...) creates a new PCA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCA_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCA_GUI

% Last Modified by GUIDE v2.5 08-Nov-2015 16:56:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCA_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PCA_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PCA_GUI is made visible.
function PCA_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCA_GUI (see VARARGIN)

% Choose default command line output for PCA_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCA_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PCA_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function nDataPairs_Callback(hObject, eventdata, handles)
% hObject    handle to nDataPairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nDataPairs as text
%        str2double(get(hObject,'String')) returns contents of nDataPairs as a double


% --- Executes during object creation, after setting all properties.
function nDataPairs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nDataPairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetCorr_Callback(hObject, eventdata, handles)
% hObject    handle to targetCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetCorr as text
%        str2double(get(hObject,'String')) returns contents of targetCorr as a double


% --- Executes during object creation, after setting all properties.
function targetCorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DC_x_Callback(hObject, eventdata, handles)
% hObject    handle to DC_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DC_x as text
%        str2double(get(hObject,'String')) returns contents of DC_x as a double


% --- Executes during object creation, after setting all properties.
function DC_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DC_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DC_y_Callback(hObject, eventdata, handles)
% hObject    handle to DC_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DC_y as text
%        str2double(get(hObject,'String')) returns contents of DC_y as a double


% --- Executes during object creation, after setting all properties.
function DC_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DC_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_data.
function update_data_Callback(hObject, eventdata, handles)
%% This is the most important function of this GUI.
% It generates the data, computes a linear fit, computes PCA, and plots everything

% -----------------------------------------------------------------------
% STEP 1: Gather parameters from other fields via 'handles'
% (the text field returns a string, and we convert to double)
nPairs = str2double( get(handles.nDataPairs,'String') ); 
datCor = str2double( get(handles.targetCorr,'String') );
DC_x   = str2double( get(handles.DC_x,'String') );
DC_y   = str2double( get(handles.DC_y,'String') );


% -----------------------------------------------------------------------
% STEP 2: CREATE NEW RANDOM DATA
% random data with forced correlation
data = randn(nPairs,2);
data(:,2) = datCor*data(:,1) + sqrt(1-datCor^2)*data(:,2);

% add X and Y offsets
data(:,1) = data(:,1) + DC_x;
data(:,2) = data(:,2) + DC_y;

% get min/max XY values for setting plot limits
plotLims = max(abs([min(data(:)) max(data(:))]));
plotLims = [-plotLims plotLims]*1.2;

% -----------------------------------------------------------------------
% STEP 3: compute linear fit and predicted values
A = [ ones(nPairs,1) data(:,1) ];
x = (A'*A)\(A'*data(:,2));
yHat = x(1)*A(:,1) + x(2)*A(:,2);


% -----------------------------------------------------------------------
% STEP 4: plot in original data axis
cla(handles.orig_axis)
plot(handles.orig_axis,data(:,1),data(:,2),'bo','markerfacecolor','k','markersize',8)
hold(handles.orig_axis,'on')
plot(handles.orig_axis,data(:,1),yHat,'r','linew',3)
set(handles.orig_axis,'xlim',plotLims,'ylim',plotLims)
grid(handles.orig_axis,'minor')
xlabel(handles.orig_axis,'X-axis'), ylabel(handles.orig_axis,'Y-axis')
legend(handles.orig_axis,{[ 'Data with r=' num2str(datCor) ];'best-fit line'})

% thick black lines on axes
plot(handles.orig_axis,[0 0],plotLims,'k','linew',2)
plot(handles.orig_axis,plotLims,[0 0],'k','linew',2)

% -----------------------------------------------------------------------
% STEP 5: run PCA
% optional subtract mean
if get(handles.meanSub,'Value')==1
    data = bsxfun(@minus,data,mean(data,1));
end
covmat = (data'*data) / nPairs;
[evecs,evals] = eig(covmat);


% -----------------------------------------------------------------------
% STEP 6: plot data in PC space 
pc1 = data*evecs(1,:)';
pc2 = data*evecs(2,:)';

cla(handles.pc_axis)
plot(handles.pc_axis,pc2,pc1,'ro','markerfacecolor','k','markersize',10)
hold on
set(handles.pc_axis,'xlim',plotLims,'ylim',plotLims)
grid(handles.pc_axis,'minor')
xlabel(handles.pc_axis,'PC1'), ylabel(handles.pc_axis,'PC2')

% thick black lines on axes
plot(handles.pc_axis,[0 0],plotLims,'k','linew',2)
plot(handles.pc_axis,plotLims,[0 0],'k','linew',2)

% --- Executes on button press in save_orig.
function save_orig_Callback(hObject, eventdata, handles)
% hObject    handle to save_orig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_orig


% --- Executes on button press in save_pc.
function save_pc_Callback(hObject, eventdata, handles)
% hObject    handle to save_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_pc


% --- Executes on button press in save2file.
function save2file_Callback(hObject, eventdata, handles)
% hObject    handle to save2file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save2matlab.
function save2matlab_Callback(hObject, eventdata, handles)
% hObject    handle to save2matlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in meanSub.
function meanSub_Callback(hObject, eventdata, handles)
% hObject    handle to meanSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meanSub
