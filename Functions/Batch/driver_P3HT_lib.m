addpath(genpath(pwd))

% Load table with processing info
ss2=import_processing('~/Google Drive/My AFM/All Good Images/All Sizes/sum_stats_revised_2.csv');
proc_table=ss2(:,{'ImageName','Process','Depo','ImageSize'});

Results_Dir = '~/Google Drive/My AFM/All Good Images/All Sizes/Unlabeled/';
BlueVec_Dir = '~/Google Drive/My AFM/All Good Images/All Sizes/BlueVec/';

dir_5_adapt = '~/Google Drive/My AFM/All Good Images/5 um/Adaptive/';
dir_5_global = '~/Google Drive/My AFM/All Good Images/5 um/Global/';
dir_4 = '~/Google Drive/My AFM/All Good Images/4 um/';
dir_7 = '~/Google Drive/My AFM/All Good Images/7 um/';

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
    'maxCurv', 7, ...
    'minFibLen', 100 ...
    );

% We would like to analyze all of the images, and save the results with
% anonymyzed names

for i = 1:height(proc_table)
    disp('_______')
    disp([num2str(i) ' of ' num2str(height(proc_table))])
    imName = proc_table(i,'ImageName');
    imCode = num2str(i);
    newFile = [Results_Dir, imCode, '.tif'];
    blueFile = [BlueVec_Dir, imCode, '_vec.tif'];
    
    % Try each directory, determine image width and whether or not to use
    % global threshold based on directory it's found in
    found = 0;
    try
        imread([dir_5_adapt,proc_table{i,'ImageName'}{1}]);
        imFile = [dir_5_adapt,proc_table{i,'ImageName'}{1}];
        Params.nmWid = 5000;
        Params.threshMethod = 1;
        found = found+1;
    catch
    end
    
    try
        imread([dir_5_global,proc_table{i,'ImageName'}{1}]);
        imFile = [dir_5_global,proc_table{i,'ImageName'}{1}];
        Params.nmWid = 5000;
        Params.threshMethod = 2;
        found = found+1;
    catch
    end
    
    try
        imread([dir_7,proc_table{i,'ImageName'}{1}]);
        imFile = [dir_7,proc_table{i,'ImageName'}{1}];
        Params.nmWid = 7000;
        Params.threshMethod = 1;
        found = found+1;
    catch
    end
    
    try
        imread([dir_4,proc_table{i,'ImageName'}{1}]);
        imFile = [dir_4,proc_table{i,'ImageName'}{1}];
        Params.nmWid = 4000;
        Params.threshMethod = 1;
        found = found+1;
    catch
    end
    
    if found ~= 1
        disp(imName)
        disp('was found in 0 or 2+ directories')
        continue
    end
    
    ims = get_ims_nogui(imFile,Params);
    [~, Fim] = FiberVecPlot_stitch(struct('ims',ims),0);
    imwrite(ims.img,newFile);
    imwrite(Fim,blueFile);
    
    proc_table{i,'ImNum'} = i;
    proc_table{i,'NewFile'} = {newFile};
    proc_table{i,'Sfull'} = ims.op2d.Sfull;
    proc_table{i,'CorrLen'} = ims.op2d.decayLen;
    proc_table{i,'AvgOrient'} = ims.ODist.director;
    proc_table{i,'LengthDensity'} = ims.fibLengthDensity;
    proc_table{i,'MeanLength'} = mean(ims.FLD);
    proc_table{i,'MeanWidth'} = mean(ims.FWD);
    proc_table{i,'a'} = ims.op2d_fa.a;
    proc_table{i,'b'} = ims.op2d_fa.b;
    proc_table{i,'lambda'} = ims.op2d_fa.lambda;
    [~,place] = min(abs(ims.op2d.xdata-5000));
    proc_table{i,'S2D5um'} = ims.op2d.S_im(place);
    proc_table{i,'a_new'} = ims.op2d.a;
    proc_table{i,'lp'} = ims.lp;
    
    close all
end

for i = 1:length(proc_table.Process);
proc_table.Process{i}=strrep(proc_table.Process{i},'2MP','PS');
end
csv_path = '/Users/Imperssonator/Google Drive/My AFM/All Good Images/All Sizes/sum_stats_revised_model.csv';
writetable(proc_table,csv_path);

