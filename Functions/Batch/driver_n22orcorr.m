im_dir = '~/CC/N2200/data/AFM/lo hi spin';
res_dir = fullfile(im_dir,'res');
fig_dir = fullfile(im_dir,'img');
nmWid = 10000;

if ~exist(res_dir, 'dir')
  mkdir(res_dir);
end
if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% Build directory of images
ad=pwd;
cd(im_dir)
dd = dir('*.jpg');
cd(ad);

% Define parameters as you would in the GUI
Params = struct( ...
    'nmWid', nmWid, ...
    'gaussnm', 10, ...
    'rhonm', 60, ...
    'T', 5, ...
    'thnm', 80, ...
    'threshMethod', 1, ...
    'globalThresh', 0.45, ...
    'noisenm', 3000, ...
    'maxBranchSizenm', 80, ...
    'fiberStep_nm', 60, ...
    'stitchGap', 80, ...
    'maxCurv', 5, ...
    'minFibLen', 0 ...
    );

% Analyze images and save results

Res = zeros(length(dd),5);

for i = 1:length(dd)
    disp('_______')
    disp([num2str(i) ' of ' num2str(length(dd))])
    imName = dd(i).name;
    imFile = fullfile(dd(i).folder,dd(i).name);
    

    blue_file = fullfile(fig_dir, [imName(1:end-4), '_vec.png']);           % File name for blue vector plot
    acm_file = fullfile(fig_dir, [imName(1:end-4), '_acm.png']);            % File name for angle color map
    s2d_file = fullfile(res_dir, [imName(1:end-4), '_s2d.png']);            % Film name for S2D plot
    orcorr_file = fullfile(res_dir, [imName(1:end-4), '_ocf2d.png']);       % File name for OrCorr plot
    odist_file = fullfile(res_dir, [imName(1:end-4), '_OD.png']);           % File name for O Dist
    mat_file = fullfile(im_dir, [imName(1:end-4), '_ims']);                % Save the ims structure

    % Process the image
    ims = get_ims_nogui(imFile,Params);
    ims = orcorr2d(ims);
    save(mat_file, ims);
    
    % Save blue vectors and plot images
    [hblue, blue_im] = FiberVecPlot_stitch(struct('ims',ims),0);
    imwrite(blue_im,blue_file);
    
    [hacm, acm_im] = FiberVec_ACM(ims,0);
    imwrite(acm_im,acm_file);
    
    [hs2d, s2d_im] = plotS2D(ims,0);
    imwrite(s2d_im,s2d_file);
    
    [hod, od_im] = ODist_plot(ims,0);
    imwrite(od_im,odist_file);
    
    [~, orcorr_im] = plotOrCorr2D(ims,0);
    imwrite(orcorr_im,orcorr_file);
    
    Res(i,:) = [ims.op2d.Sfull, ...
                ims.op2d.decayLen, ...
                ims.ODist.director, ...
                ims.fibLengthDensity, ...
                ims.op2d.a];
    
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
end

close all

dd=rmfield(dd,'folder');
csv_path = fullfile(res_dir, 'gtfiber_results.csv');
struct2csv(dd,csv_path);