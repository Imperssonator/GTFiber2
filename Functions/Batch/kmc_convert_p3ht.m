im_dir = 'simulation/p3ht';
fig_dir = fullfile(im_dir,'img');
nmWid = 5000;

if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% Build directory of images
ad=pwd;
cd(im_dir)
dd = dir('*.tif');
cd(ad);

% Define parameters as you would in the GUI
Params = struct( ...
    'nmWid', nmWid, ...
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
    'maxCurv', 7, ...
    'minFibLen', 0 ...
    );

% Analyze images and save results

% Res = zeros(length(dd),6);

for i = 1:length(dd)
    disp('_______')
    disp([num2str(i) ' of ' num2str(length(dd))])
    imName = dd(i).name;
    imFile = fullfile(dd(i).folder,dd(i).name);
    

    bw_file = fullfile(fig_dir, [imName(1:end-4), '_bw.png']);           % File name for binarized image
    bw_array_file = fullfile(fig_dir, [imName(1:end-4), '_bw_array.txt']);   % File name for binarized array
    acm_file = fullfile(fig_dir, [imName(1:end-4), '_acm.png']);          % File name for angle color map
    acm_array_file = fullfile(fig_dir, [imName(1:end-4), '_acm_array.txt']);   % File name for binarized array
    mat_file = fullfile(im_dir, [imName(1:end-4), '_ims.mat']);           % Save the ims structure

    % Process the image
    ims = get_ims_nogui_nostitch(imFile,Params);
    
%     ims = orcorr2d(ims);
    save(mat_file, 'ims');
    
    % Save binarized image as image and ascii
    imwrite(ims.CEDclean,bw_file);
    dlmwrite(bw_array_file,ims.CEDclean);
    
    acm_im = AngleColorImage(ims.AngMap,ims.CEDclean);
    imwrite(acm_im,acm_file);
    dlmwrite(acm_array_file,ims.AngMap);    % AngMap is angle off the image horizontal
    
%     [hs2d, s2d_im] = plotS2D(ims,1);
%     [hs2d, s2d_im] = plotS2D_fixed(ims,1);
%     [hod, od_im] = ODist_plot(ims,1);
% %     [~, orcorr_im] = plotOrCorr2D(ims,1);
%     
%     Res(i,:) = [ims.op2d.Sfull, ...
%                 ims.op2d.decayLen, ...
%                 ims.ODist.director, ...
%                 ims.fibLengthDensity, ...
%                 ims.op2d.a, ...
%                 ims.op2d_fixed.S2D];
%     
end

% for i = 1:length(dd)
%     dd(i).Sfull = Res(i,1);
%     dd(i).decayLen = Res(i,2);
%     dd(i).avgOrient = Res(i,3);
%     dd(i).lengthDensity = Res(i,4);
%     dd(i).a = Res(i,5);
%     dd(i).S2D = Res(i,6);
% end

close all

% dd=rmfield(dd,'folder');
% csv_path = fullfile(im_dir, 'gtfiber_results.csv');
% struct2csv(dd,csv_path);