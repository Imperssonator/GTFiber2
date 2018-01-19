im_dir = 'F:\GitHub\basic-afm\data\171228 NDI log speed\phase';
results_dir = 'F:\GitHub\basic-afm\data\171228 NDI log speed\gtfiber\res';
blue_vec_dir = 'F:\GitHub\basic-afm\data\171228 NDI log speed\gtfiber\blue';

% Build directory of images
ad=pwd;
cd(im_dir)
dd = dir('*.png');
cd(ad);

% Define parameters as you would in the GUI
Params = struct( ...
    'nmWid', 5000, ...
    'gaussnm', 5, ...
    'rhonm', 15, ...
    'T', 5, ...
    'thnm', 40, ...
    'threshMethod', 1, ...
    'globalThresh', 0.45, ...
    'noisenm', 1500, ...
    'maxBranchSizenm', 60, ...
    'fiberStep_nm', 30, ...
    'stitchGap', 60, ...
    'maxCurv', 8, ...
    'minFibLen', 0 ...
    );

% Analyze images and save results

Res = zeros(length(dd),6);

parfor i = 1:length(dd)
    disp('_______')
    disp([num2str(i) ' of ' num2str(length(dd))])
    imName = dd(i).name;
    imFile = fullfile(dd(i).folder,dd(i).name);
    
    % File name for blue vector plot
    blue_file = fullfile(blue_vec_dir, [imName(1:end-4), '_vec.png']);
    % Film name for S2D plot
    s2d_file = fullfile(results_dir, [imName(1:end-4), '_s2d.png']);
    % File name for OrCorr plot
    orcorr_file = fullfile(results_dir, [imName(1:end-4), '_orcorr.png']);

    % Process the image
    ims = get_ims_nogui(imFile,Params);
    
    % Save blue vectors and plot images
    [~, blue_im] = FiberVecPlot_stitch(struct('ims',ims),0);
    imwrite(blue_im,blue_file);
    
    [~, s2d_im] = plotS2D(ims,0);
    imwrite(s2d_im,s2d_file);
    
    [~, orcorr_im] = plotOrCorr(ims,0);
    imwrite(orcorr_im,orcorr_file);
    
    Res(i,:) = [ims.op2d.Sfull, ...
                ims.op2d.decayLen, ...
                ims.ODist.director, ...
                ims.fibLengthDensity, ...
                ims.op2d.a, ...
                ims.OrCorr.corrLen];
    
%     dd(i).Sfull = ims.op2d.Sfull;
%     dd(i).decayLen = ims.op2d.decayLen;
%     dd(i).avgOrient = ims.ODist.director;
%     dd(i).lengthDensity = ims.fibLengthDensity;
%     dd(i).a = ims.op2d.a;
%     dd(i).corrLen = ims.OrCorr.corrLen;
    
end

for i = 1:length(dd)
    dd(i).Sfull = Res(i,1);
    dd(i).decayLen = Res(i,2);
    dd(i).avgOrient = Res(i,3);
    dd(i).lengthDensity = Res(i,4);
    dd(i).a = Res(i,5);
    dd(i).corrLen = Res(i,6);
end

close all

dd=rmfield(dd,'folder');
csv_path = fullfile(results_dir, 'gtfiber_results.csv');
struct2csv(dd,csv_path);