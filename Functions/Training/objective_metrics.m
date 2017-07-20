function [out, results] = objective_metrics(Params,imFile,nmWid,fa_path)

% This objective function calculates the Euclidean distance of the
% structural metrics from GTFiber from the structural metrics from a manual
% tracing in FiberApp

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

f = figure;
auto_handles.img_axes = axes();
f2 = figure;
auto_handles.fiber_axes = axes();

auto_handles = main_filter(auto_handles);
ims = auto_handles.ims;

ims = StitchFibers2(ims);

[out,metrics] = metric_norm(fa_path,ims);
disp('score')
disp(out)

results = ...
struct(...
'Params',Params,...
'Fibers',ims.Fibers,...
'fibLengthDensity',ims.fibLengthDensity,...
'op2d',ims.op2d,...
'ODist',ims.ODist,...
'FLD',ims.FLD,...
'score',metrics);

close(f)
close(f2)

end

function [out,gt_metrics] = metric_norm(fa_path,ims)

% Get GT Measures
gt_metrics = [ ...
    ims.op2d.Sfull;
    ims.op2d.decayLen;
    mean(ims.FLD);
    median(ims.FLD);
    ims.fibLengthDensity];

% Get FA measures
fld_fa=FiberLengthsFA(fa_path);
load(fa_path)
[op_fa,fibLenDens]=op2d_FA(imageData);
fa_metrics = [ ...
    op_fa.Sfull;
    op_fa.decayLen;
    mean(fld_fa);
    median(fld_fa);
    fibLenDens];

% Score is the Euclidean distance of each metric away from its ground truth
% normalized by the value of the ground truth. (i.e. the Euclidean distance
% in percentage-from-truth space)
out = norm((gt_metrics-fa_metrics)./fa_metrics);

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
settings.globalThresh = 0; %str2num(get(handles.globalThresh,'String'));
settings.threshMethod = 1; %get(handles.threshMethod,'Value');

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