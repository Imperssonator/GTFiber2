% im_dir = '/Users/nils/CC/N2200/data/AFM/180726 N2200 Subh/png/Phase';
im_dir = 'Z:\AFM\181109 DPP Speed Tests\phase';
fig_dir = fullfile(im_dir,'img');
nmWid = 10000;

if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% Build directory of images
ad=pwd;
cd(im_dir)
dd = dir('*.png');
cd(ad);

dd=dd(10:11);

% Define parameters as you would in the GUI
Params = struct( ...
    'nmWid', nmWid, ...
    'gaussnm', 10, ...
    'rhonm', 30, ...
    'T', 3, ...
    'thnm', 200, ...
    'threshMethod', 1, ...
    'globalThresh', 0.45, ...
    'noisenm', 3000, ...
    'maxBranchSizenm', 120, ...
    'fiberStep_nm', 60, ...
    'stitchGap', 80, ...
    'maxCurv', 5, ...
    'minFibLen', 0 ...
    );

% Analyze images and save results

Res = zeros(length(dd),3);

for i = 1:length(dd)
    disp('_______')
    disp([num2str(i) ' of ' num2str(length(dd))])
    
    imName = dd(i).name;
    imFile = fullfile(dd(i).folder,dd(i).name);
    
    acm_file = fullfile(fig_dir, [imName(1:end-4), '_acm.png']);            % File name for angle color map
    mat_file = fullfile(im_dir, [imName(1:end-4), '_ims.mat']);                % Save the ims structure

    % Process the image
    ims = get_ims_nogui_nostitch(imFile,Params);
    ims = ODist_noVec(ims);
    ims = orcorr2d_noVec(ims);
    
%     ims = orcorr2d(ims);
    save(mat_file, 'ims');
    
    % Save blue vectors and plot images
    acm_im = AngleColorImage(ims.AngMap, ims.CEDclean);
    imwrite(acm_im,acm_file);
    
    [hod, od_im] = ODist_plot(ims,1);
    [~, orcorr_im] = plotOrCorr2D(ims,1);
    
    Res(i,:) = [ims.ODist.S, ...
                ims.OrCorr2D.corrLen_x, ...
                ims.OrCorr2D.corrLen_y];
    
%     dd(i).Sfull = ims.op2d.Sfull;
%     dd(i).decayLen = ims.op2d.decayLen;
%     dd(i).avgOrient = ims.ODist.director;
%     dd(i).lengthDensity = ims.fibLengthDensity;
%     dd(i).a = ims.op2d.a;
%     dd(i).corrLen = ims.OrCorr.corrLen;
    
end

for i = 1:length(dd)
    dd(i).S2D = Res(i,1);
    dd(i).corrLen_x = Res(i,2);
    dd(i).corrLen_y = Res(i,3);
end

close all

dd=rmfield(dd,'folder');
csv_path = fullfile(im_dir, 'gtfiber_results_2.csv');
struct2csv(dd,csv_path);