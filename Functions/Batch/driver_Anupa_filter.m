im_dir = '~/CC/anupa/filt';
results_dir = fullfile(im_dir,'res');
blue_vec_dir = fullfile(im_dir,'img');

if ~exist(results_dir, 'dir')
  mkdir(results_dir);
end
if ~exist(blue_vec_dir, 'dir')
  mkdir(blue_vec_dir);
end

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
    'maxBranchSizenm', 20, ...
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
%     orcorr_file = fullfile(results_dir, [imName(1:end-4), '_orcorr.png']);
    % File name for labeled nodes
    node_file = fullfile(blue_vec_dir, [imName(1:end-4), '_nodes.png']);
    % File name for O Dist
    odist_file = fullfile(results_dir, [imName(1:end-4), '_OD.png']);

    % Process the image
    ims = get_ims_nogui(imFile,Params);
    
    % Save blue vectors and plot images
%     [~, blue_im] = FiberVecPlot_stitch(struct('ims',ims),0);
    imwrite(ims.CEDclean,blue_file);
    
    [~, s2d_im] = plotS2D(ims,0);
    imwrite(s2d_im,s2d_file);
    
    [bpd, bp_im] = node_analysis(ims)
    imwrite(bp_im,node_file);
    
    [~, od_im] = ODist_plot(ims,0);
    imwrite(od_im,odist_file);
    
%     [~, orcorr_im] = plotOrCorr(ims,0);
%     imwrite(orcorr_im,orcorr_file);
    
    Res(i,:) = [ims.op2d.Sfull, ...
                ims.op2d.decayLen, ...
                ims.ODist.director, ...
                ims.fibLengthDensity, ...
                ims.op2d.a, ...
                bpd];
    
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
    dd(i).numNodes = Res(i,6);
end

close all

dd=rmfield(dd,'folder');
csv_path = fullfile(results_dir, 'gtfiber_results.csv');
struct2csv(dd,csv_path);