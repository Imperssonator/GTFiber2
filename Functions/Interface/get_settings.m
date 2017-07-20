function settings = get_settings(handles)

% Get dimensions
settings.nmWid = str2num(get(handles.nmWid,'String'));

% Get filter settings
settings.invert = get(handles.invertColor,'Value');
settings.thnm = str2num(get(handles.tophatSize,'String'));
settings.noisenm = str2num(get(handles.noiseArea,'String'));
settings.maxBranchSizenm = str2num(get(handles.maxBranchSize,'String'));
% settings.maxStubLennm = str2num(get(handles.maxStubLen,'String'));
settings.globalThresh = str2num(get(handles.globalThresh,'String'));
settings.threshMethod = get(handles.threshMethod,'Value');

settings.figSave = get(handles.saveFigs,'Value');

% Build the Coherence Filter options structure
Options = struct();
Options.Scheme = 'I';
settings.gaussnm = str2num(get(handles.gauss,'String'));
settings.rhonm = str2num(get(handles.rho,'String'));
% Options.sigma = gausspix;
% Options.rho = rhopix;
Options.T = str2num(get(handles.difftime,'String'));
Options.dt = 0.15;
% Options.eigenmode = 5;
Options.eigenmode = 0;
Options.C = 1E-10;

settings.Options = Options;


% OP2D calculation settings
% settings.gridStepnm = str2num(get(handles.gridStep,'String'));
% settings.frameStepnmWide = str2num(get(handles.frameStep,'String'));

% Fiber Width calc settings
% settings.fibWidSamps = str2num(get(handles.fibWidSamps,'String'));
settings.fibWidSamps2 = 15;


% Fiber Fitting Settings
settings.fiberStep_nm = str2num(get(handles.fiberStep,'String'));
% settings.maxAngleDeg = str2num(get(handles.maxAngleDeg,'String'));
settings.maxCurv = str2num(get(handles.maxCurv,'String')) / 1000;
% settings.curvLen = str2num(get(handles.curvLen,'String'));
% settings.minWidthnm = str2num(get(handles.minWidth,'String'));
% settings.maxWidthnm = str2num(get(handles.maxWidth,'String'));
settings.stitchGap = str2num(get(handles.stitchGap,'String'));
settings.minFibLen = str2num(get(handles.minFibLen,'String'));


% Gif Export Settings
settings.initDelay = 1;
settings.CEDStepDelay = 0.1;
settings.CEDFinalDelay = 0.8;
settings.skelDelay = 4;
settings.plotDelay = 0.5;
settings.plotFinal = 3;

end