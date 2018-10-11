function ims = get_ims_nogui_nostitch(imFile,Params)

% imFile is the path to an image file
% Params is a struct with fields:
%
% nmWid (nm)
% gaussnm (nm)
% rhonm (nm)
% T (nm)
% thnm (nm)
% threshMethod (1=adaptive, 2=global)
% globalThresh (0-1)
% noisenm (nm^2)
% maxBranchSizenm (nm)
% fiberStep_nm (nm)
% stitchGap (nm)
% maxCurv (1/um)
% minFibLen (nm)

addpath(genpath(pwd))

% Initialize image data
ims = initImgData(imFile);

% Build settings from the provided "Params"
ims.settings = get_settings_nogui(Params);
ims = pix_settings(ims);

% Make a dummy handles structure
auto_handles.ims = ims;

% Run filter and stitch fibers
auto_handles = main_filter(auto_handles);
ims = auto_handles.ims;

end


function settings = get_settings_nogui(Params)

% Get dimensions
settings.nmWid = Params.nmWid; %str2num(get(handles.nmWid,'String'));
settings.invert = 0; %get(handles.invertColor,'Value');

% Build the Coherence Filter options structure
Options = struct();
Options.Scheme = 'I';
settings.gaussnm = Params.gaussnm; %str2num(get(handles.gauss,'String'));
settings.rhonm = Params.rhonm; %str2num(get(handles.rho,'String'));
Options.T = Params.T; %str2num(get(handles.difftime,'String'));
Options.dt = 0.15;
% Options.eigenmode = 5;
Options.eigenmode = 0;
Options.C = 1E-10;

settings.Options = Options;

% Top Hat Filter
settings.thnm = Params.thnm; %str2num(get(handles.tophatSize,'String'));

% Thresholding
settings.globalThresh = Params.globalThresh; %str2num(get(handles.globalThresh,'String'));
settings.threshMethod = Params.threshMethod; %get(handles.threshMethod,'Value');
settings.noisenm = Params.noisenm; %str2num(get(handles.noiseArea,'String'));

% Skeletonization
settings.maxBranchSizenm = Params.maxBranchSizenm; %str2num(get(handles.maxBranchSize,'String'));

% Stitching
settings.fiberStep_nm = Params.fiberStep_nm;
settings.stitchGap = Params.stitchGap;
settings.maxCurv = Params.maxCurv / 1000;   %!!!
settings.minFibLen = Params.minFibLen;

% Batch Processing
settings.figSave = 0; %get(handles.saveFigs,'Value');
settings.figSwitch = 0;

% Fiber Width calc settings
settings.fibWidSamps2 = 15;

% Gif Export Settings
settings.initDelay = 1;
settings.CEDStepDelay = 0.1;
settings.CEDFinalDelay = 0.8;
settings.skelDelay = 4;
settings.plotDelay = 0.5;
settings.plotFinal = 3;

end