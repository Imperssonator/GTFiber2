% Copyright (C) 2016 Nils Persson
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function varargout = GTFiber(varargin)
% GTFIBER MATLAB code for GTFiber.fig
%      GTFIBER, by itself, creates a new GTFIBER or raises the existing
%      singleton*.
%
%      H = GTFIBER returns the handle to a new GTFIBER or the handle to
%      the existing singleton*.
%
%      GTFIBER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GTFIBER.M with the given input arguments.
%
%      GTFIBER('Property','Value',...) creates a new GTFIBER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GTFiber_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GTFiber_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GTFiber

% Last Modified by GUIDE v2.5 16-Jun-2017 18:55:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GTFiber_OpeningFcn, ...
                   'gui_OutputFcn',  @GTFiber_OutputFcn, ...
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


function GTFiber_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for GTFiber
handles.output = hObject;
% addpath(genpath('./'))

% Update handles structure
guidata(hObject, handles);


function varargout = GTFiber_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%__________________________________________________________________________


function Main_Callback(hObject, eventdata, handles)


function Load_Callback(hObject, eventdata, handles)

[filename, folderpath] = uigetfile({'*.jpg;*.jpeg;*.tif;*.tiff;*.png;*.gif;*.bmp','All Image Files'});
if isequal(filename, 0); return; end % Cancel button pressed

% Prompt user for the image width
prompt = {'Enter image width in nanometers, with no commas (ex. 5000):'};
dlg_title = 'Image Scale';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
set(handles.nmWid,'String',answer{1})
nmWid_Callback(hObject, eventdata, handles);

% Initialize the internal image data structure, "ims"
imfile = [folderpath, filename];
handles.ims = initImgData(imfile);
set(handles.fileNameBox,'String',handles.ims.imName);

% Initialize the figure window and don't let the user close it
handles = imshowGT(handles.ims.img,handles,'img_axes');

guidata(hObject, handles);


function Coherence_Filter_Callback(hObject, eventdata, handles)

if ~isfield(handles,'ims')
    noload = errordlg('Go to File>Load Image to load an image before filtering.');
    return
end

% Get Settings
handles.ims.settings = get_settings(handles);
handles.ims = pix_settings(handles.ims);

% Run Filter Regime
handles = main_filter(handles);

guidata(hObject, handles);


function runStitch_Callback(hObject, eventdata, handles)

if ~isfield(handles.ims,'skelTrim')
    noload = errordlg('Go to File>Load Image to load an image, then Run Filter.');
    return
end

% Get Settings
handles.ims.settings = get_settings(handles);
handles.ims = pix_settings(handles.ims);

% Stitch fiber segments and calculate length
handles.ims = StitchFibers2(handles.ims);
handles = FiberVecPlot_stitch(handles);

guidata(hObject, handles);


function AngMap_Callback(hObject, eventdata, handles)

if ~isfield(handles,'ims')
    noload = errordlg('Go to File>Load Image to load an image before filtering.');
    return
end

if ~isfield(handles.ims,'AngMap')
    nofilt = errordlg('"Run Filter" must be executed before results can be displayed');
    return
end

FiberVec_ACM(handles.ims);


function op2d_Callback(hObject, eventdata, handles)

if ~isfield(handles,'ims')
    noload = errordlg('Go to File>Load Image to load an image before filtering.');
    return
end

if ~isfield(handles.ims,'AngMap')
    nofilt = errordlg('"Run Filter" must be executed before results can be displayed');
    return
end

if ~isfield(handles.ims,'Fibers')
    nofilt = errordlg('"Stitch Fibers" must be executed before results can be displayed');
    return
end

plotS2D(handles.ims,0);
ODist_plot(handles.ims,0);

guidata(hObject, handles);


function GetFiberLength_Callback(hObject, eventdata, handles)

if ~isfield(handles,'ims')
    noload = errordlg('Go to File>Load Image to load an image before filtering.');
    return
end

if ~isfield(handles.ims,'Fibers')
    nofilt = errordlg('"Stitch Fibers" must be executed before results can be displayed');
    return
end

FLD_hist(handles.ims);
FWD_hist(handles.ims);

guidata(hObject, handles);


% --- Executes on button press in runDir.
function runDir_Callback(hObject, eventdata, handles)

% Solicit the folder to run
folderPath = uigetdir;

if isequal(folderPath, 0)
    return
end

if ispc
    separator = '\';
else
    separator = '/';
end

folderPath = [folderPath, separator];

% Get name for results file
prompt = {'Save results with file name (no extension necessary):'};
dlg_title = 'Save File Name';
num_lines = 1;
fileName = inputdlg(prompt,dlg_title,num_lines);
saveFilePath = [folderPath, fileName{1}, '.csv'];

run_directory(handles,folderPath,saveFilePath);


% function chainStacker_Callback(hObject, eventdata, handles)
% 
% % Solicit the folder to run
% folderPath = uigetdir;
% 
% if isequal(folderPath, 0)
%     return
% end
% 
% if ispc
%     separator = '\';
% else
%     separator = '/';
% end
% 
% folderPath = [folderPath, separator];
% 
% % Get name for results file
% prompt = {'Save results with file name (no extension necessary):'};
% dlg_title = 'Save File Name';
% num_lines = 1;
% fileName = inputdlg(prompt,dlg_title,num_lines);
% saveFilePath = [folderPath, fileName{1}, '.csv'];
% 
% % Solicit range of molecular weight and PDI for simulations
% prompt = {'List of Mn to use', 'List of PDI to use', 'Number of simulations per point'};
% dlg_title = 'MW Parameters';
% num_lines = 3;
% Mn_input = inputdlg(prompt,dlg_title,num_lines);
% 
% Mn_vec = str2num(Mn_input{1});
% PDI_vec = str2num(Mn_input{2});
% nruns = str2num(Mn_input{3});
% 
% run_directory_CS(handles,folderPath,saveFilePath,Mn_vec,PDI_vec,nruns);



%__________________________________________________________________________
% Fields for image processing settings


function gauss_Callback(hObject, eventdata, handles)

function gauss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    
end


function rho_Callback(hObject, eventdata, handles)

function rho_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function difftime_Callback(hObject, eventdata, handles)

function difftime_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noiseArea_Callback(hObject, eventdata, handles)

function noiseArea_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minFibLen_Callback(hObject, eventdata, handles)

function minFibLen_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tophatSize_Callback(hObject, eventdata, handles)

function tophatSize_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nmWid_Callback(hObject, eventdata, handles)

nmWid = str2num(get(handles.nmWid,'String'));

if get(handles.scaleParams,'Value')
    if ~isempty(nmWid)
        set(handles.gauss,'String',num2str(nmWid*5/5000));
        set(handles.rho,'String',num2str(nmWid*15/5000));
        set(handles.tophatSize,'String',num2str(nmWid*40/5000));
        set(handles.noiseArea,'String',num2str(nmWid*1500/5000));
        set(handles.maxBranchSize,'String',num2str(nmWid*60/5000));
        set(handles.stitchGap,'String',num2str(nmWid*60/5000));
    end
end

guidata(hObject, handles);

function nmWid_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function threshMethod_Callback(hObject, eventdata, handles)

function threshMethod_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function globalThresh_Callback(hObject, eventdata, handles)

function globalThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxBranchSize_Callback(hObject, eventdata, handles)

function maxBranchSize_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function saveFigs_Callback(hObject, eventdata, handles)


function mainFig_CreateFcn(hObject, eventdata, handles)


function widthText_CreateFcn(hObject, eventdata, handles)


function fibWidSamps_Callback(hObject, eventdata, handles)

function fibWidSamps_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%__________________________________________________________________________


% function Make_Gif_Callback(hObject, eventdata, handles)
% if ~isfield(handles,'ims')
%     noload = errordlg('Go to File>Load Image to load an image before filtering.');
%     return
% end
% 
% settings = get_settings(handles);
% [settings, ims] = pix_settings(settings,handles.ims);
% handles.ims = ims;
% gif_filter(handles.ims,settings);
% 
% settings.figSwitch = 1; % Gotta turn on figSwitch to make the figure
% settings.figSave = 0;   % No need to save
% gif_op2d_am(handles.ims,settings);
% 
% guidata(hObject, handles);


function invertColor_Callback(hObject, eventdata, handles)

switch get(handles.invertColor,'Value')
    case 1
        handles=imshowGT(imcomplement(handles.ims.gray),handles,'img_axes');
    case 0
        handles=imshowGT(handles.ims.gray,handles,'img_axes');
end

guidata(hObject, handles);

%__________________________________________________________________________
% Buttons for switching what figure is displayed in the image processing
% preview window

function showCED_ButtonDownFcn(hObject, eventdata, handles)


function showTopHat_ButtonDownFcn(hObject, eventdata, handles)


function showThresh_ButtonDownFcn(hObject, eventdata, handles)


function showClean_ButtonDownFcn(hObject, eventdata, handles)


function showSkel_ButtonDownFcn(hObject, eventdata, handles)


function showSkelTrim_ButtonDownFcn(hObject, eventdata, handles)


function showSegs_ButtonDownFcn(hObject, eventdata, handles)


function showFibers_ButtonDownFcn(hObject, eventdata, handles)


function showImg_ButtonDownFcn(hObject, eventdata, handles)


function showImg_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'img')
    handles=imshowGT(handles.ims.img,handles,'img_axes');
end
guidata(hObject, handles);


function showCED_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'CEDgray')
    handles=imshowGT(handles.ims.CEDgray,handles,'img_axes');
end
guidata(hObject, handles);


function showTopHat_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'CEDtophat')
    handles=imshowGT(handles.ims.CEDtophat,handles,'img_axes');
end
guidata(hObject, handles);


function showThresh_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'CEDbw')
    handles=imshowGT(handles.ims.CEDbw,handles,'img_axes');
end
guidata(hObject, handles);


function showClean_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'CEDclean')
    handles=imshowGT(handles.ims.CEDclean,handles,'img_axes');
end
guidata(hObject, handles);


function showSkel_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'skel')
    handles=imshowGT(handles.ims.skel,handles,'img_axes');
end
guidata(hObject, handles);


function showSkelTrim_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'skelTrim')
    handles=imshowGT(handles.ims.skelTrim,handles,'img_axes');
end
guidata(hObject, handles);


function showSegs_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'fibSegs')
    handles=FiberVecPlot_segs(handles);
end
guidata(hObject, handles);


function showFibers_Callback(hObject, eventdata, handles)

if isfield(handles.ims,'Fibers')
    handles=FiberVecPlot_stitch(handles);
end
guidata(hObject, handles);

%__________________________________________________________________________
% Fiber Stitching Setting Fields

function curvLen_Callback(hObject, eventdata, handles)

function curvLen_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxCurv_Callback(hObject, eventdata, handles)

function maxCurv_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minWidth_Callback(hObject, eventdata, handles)

function minWidth_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function maxWidth_Callback(hObject, eventdata, handles)

function maxWidth_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fiberStep_Callback(hObject, eventdata, handles)

function fiberStep_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stitchGap_Callback(hObject, eventdata, handles)

function stitchGap_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function scaleParams_Callback(hObject, eventdata, handles)



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to minFibLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minFibLen as text
%        str2double(get(hObject,'String')) returns contents of minFibLen as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minFibLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
