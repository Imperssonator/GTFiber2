function settings = get_settings(handles)

% Get dimensions
settings.nmWid = str2num(get(handles.nmWid,'String'));
settings.invert = get(handles.invertColor,'Value');

% Diffusion Filter Options
Options = struct();
Options.Scheme = 'I';
settings.gaussnm = str2num(get(handles.gauss,'String'));
settings.rhonm = str2num(get(handles.rho,'String'));
Options.T = str2num(get(handles.difftime,'String'));
Options.dt = 0.15;
% Options.eigenmode = 5;
Options.eigenmode = 0;
Options.C = 1E-10;

settings.Options = Options;

% Top Hat Filter settings
settings.thnm = str2num(get(handles.tophatSize,'String'));

% Thresholding
settings.threshMethod = get(handles.threshMethod,'Value');
settings.globalThresh = str2num(get(handles.globalThresh,'String'));
settings.noisenm = str2num(get(handles.noiseArea,'String'));

% Skeletonization
settings.maxBranchSizenm = str2num(get(handles.maxBranchSize,'String'));

% Fiber Fitting Settings
settings.fiberStep_nm = str2num(get(handles.fiberStep,'String'));
settings.stitchGap = str2num(get(handles.stitchGap,'String'));
settings.maxCurv = str2num(get(handles.maxCurv,'String')) / 1000;
settings.minFibLen = str2num(get(handles.minFibLen,'String'));

% Batch Processing
settings.figSave = get(handles.saveFigs,'Value');

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