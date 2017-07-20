function results = get_results(Params,imFile,nmWid)

% Params:                                   LB      |   UB     
% 
% 1. Gaussian Smoothing (nm)                1       |   100
% 2. Orientation Smoothing (nm)             1       |   100
% 3. Diffusion Time (s)                     0.5     |   10
% 4. Top Hat Size (nm)                      1       |   100
% 5. Noise Removal (nm^2)                   100     |   5000
% 6. Fringe Removal (nm)                    0       |   300
% 7. Max Stitch Gap (nm)                    0       |   500
% 8. Max Curvature (1/pixel)                0       |   0.2
% 9. Step Length (pix)                      2       |   10

addpath(genpath(pwd))

ims = initImgData(imFile);
settings = get_settings_nogui(Params,nmWid);
settings.figSwitch = 0;
settings.figSave = 0;
settings.fullOP = 1;
[settings, ims] = pix_settings_nogui(settings,ims);
ims.settings = settings;

auto_handles.ims = ims;

auto_handles = main_filter(auto_handles);
ims = auto_handles.ims;

ims = StitchFibers2(ims);

results = ...
struct(...
'Params',Params,...
'Fibers',ims.Fibers,...
'fibLengthDensity',ims.fibLengthDensity,...
'op2d',ims.op2d,...
'ODist',ims.ODist,...
'FLD',ims.FLD);


end


function settings = get_settings_nogui(Params,nmWid)

% Get dimensions
settings.nmWid = nmWid; %str2num(get(handles.nmWid,'String'));

% Get filter settings
settings.invert = 0; %get(handles.invertColor,'Value');
settings.thnm = Params(4); %str2num(get(handles.tophatSize,'String'));
settings.noisenm = Params(5); %str2num(get(handles.noiseArea,'String'));
settings.maxBranchSizenm = Params(6); %str2num(get(handles.maxBranchSize,'String'));
settings.stitchGap = Params(7);
settings.globalThresh = 0.4; %str2num(get(handles.globalThresh,'String'));
settings.threshMethod = 2; %get(handles.threshMethod,'Value');

settings.figSave = 0; %get(handles.saveFigs,'Value');

% Build the Coherence Filter options structure
Options = struct();
Options.Scheme = 'I';
settings.gaussnm = Params(1); %str2num(get(handles.gauss,'String'));
settings.rhonm = Params(2); %str2num(get(handles.rho,'String'));
Options.T = Params(3); %str2num(get(handles.difftime,'String'));
Options.dt = 0.15;
% Options.eigenmode = 5;
Options.eigenmode = 0;
Options.C = 1E-10;

settings.Options = Options;

% Fiber Width calc settings
settings.fibWidSamps2 = 15;

% Fiber Fitting Settings
settings.fiberStep = Params(9); %ceil(str2num(get(handles.fiberStep,'String')));
settings.maxCurv = Params(8); %str2num(get(handles.maxAngleDeg,'String'));

% Gif Export Settings
settings.initDelay = 1;
settings.CEDStepDelay = 0.1;
settings.CEDFinalDelay = 0.8;
settings.skelDelay = 4;
settings.plotDelay = 0.5;
settings.plotFinal = 3;

end

function [settings, ims] = pix_settings_nogui(settings,ims)

ims.nmWid = settings.nmWid;
ims.pixWid = size(ims.img,2);
ims.nmPix = ims.nmWid/ims.pixWid;

% Get pixel values for filter options
settings.thpix = ceil(settings.thnm/ims.nmPix);
settings.noisepix = ceil(settings.noisenm/ims.nmPix^2);
settings.maxBranchSize = ceil(settings.maxBranchSizenm/ims.nmPix);
settings.Options.sigma = settings.gaussnm/ims.nmPix;
settings.Options.rho = settings.rhonm/ims.nmPix;

% Match Search Settings
settings.fiberStep_nm = settings.fiberStep * ims.nmPix;
settings.searchLat = 0.04 * ims.nmWid;
settings.searchLong = 0.04 * ims.nmWid;

end