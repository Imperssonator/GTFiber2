function xl = run_directory(handles,dirPath,savePath)

% Compile images from the directory, run the full processing/analysis stack
% specified by the current GUI parameters, save results to csv and save
% visualizations to folder with image file's name

hwaitdir = waitbar(0,'Running Directory...');

imdir = CompileImgs(dirPath);
numIms = length(imdir);

% Initialize the Cell for the csv file
xl = cell(length(imdir)+1,7);
xl{1,1} = 'Image Name';
xl{1,2} = 'Sfull fit';
xl{1,3} = 'Correlation Length (nm)';
xl{1,4} = 'Average Orientation (degrees)';
xl{1,5} = 'Fiber Length Density (1/um)';
xl{1,6} = 'Mean Fiber Length (nm)';
xl{1,7} = 'Mean Fiber Width (nm)';


for i = 1:numIms
    
    waitbar(i/numIms,hwaitdir,['Processing ', imdir(i).name]);
    
    imfilei = imdir(i).path;
%     addpath(genpath(pwd));
    
    % Initialize image data structure
    ims = initImgData(imfilei);
    ims.settings = get_settings(handles);
    ims = pix_settings(ims);
    
    if ims.settings.figSave
        mkdir(ims.imNamePath);  % make the directory to save the results figures
    end
    
    % Run the filter bank at the current settings
    handles.ims=ims;
    handles = main_filter(handles);
    ims=handles.ims;
    
    % Stitch fiber segments and calculate length
    ims = StitchFibers2(ims);
    handles.ims = ims;
    
    % Write data to csv cell
    xl{i+1,1} = imdir(i).name;
    xl{i+1,2} = ims.op2d.Sfull;
    xl{i+1,3} = ims.op2d.decayLen;
    xl{i+1,4} = ims.ODist.director;
    xl{i+1,5} = ims.fibLengthDensity;
    xl{i+1,6} = mean(ims.FLD);
    xl{i+1,7} = mean(ims.FWD);
    
    % Save figures if specified
    if ims.settings.figSave
        ODist_plot(ims,1);
        plotS2D(ims,1);
        handles = FiberVecPlot_stitch(handles,1);
        FLD_hist(ims,1);
        FWD_hist(ims,1);
    end
    
end

cell2csv(savePath, xl, ',', 1999, '.');
close(hwaitdir)

end

function out = CompileImgs(FolderPath)
disp(FolderPath)

ad = pwd;

% First compile any images from the folderpath
cd(FolderPath)

PNG = dir('*.png');
JPG = dir('*.jpg');
JPEG = dir('*.jpeg');
TIF = dir('*.tif');
TIFF = dir('*.tiff');
BMP = dir('*.bmp');
CurIms = [PNG; JPG; JPEG; TIF; TIFF; BMP]; % Generate directory structure of images in FolderPath
cd(ad)

for p = 1:length(CurIms)
    CurIms(p).path = [FolderPath, CurIms(p).name];   % prepend the folder path to the image names
end

% Remove any ghost files with the ._ prefix
c=1;
pp=1;
while pp<=length(CurIms)
    if ~strcmp(CurIms(pp).name(1:2),'._')
        out(c) = CurIms(pp);
        c = c+1;
    end
    pp = pp+1;
end
end
