function xl = runDirFLD(dirPath,settings)

% varargin{1} is: 0 if just get images, 1 if get mobilities

imdir = CompileImgs(dirPath);

% S = zeros(length(imdir),1);

DataMat = [];
auto_handles = struct();

hwaitdir = waitbar(0,'Running Directory...');
numIms = length(imdir);

for i = 1:numIms
    waitbar(i/numIms,hwaitdir,['Processing ', imdir(i).name]);
    
    imfilei = imdir(i).path;
    addpath('Functions')
    addpath('Functions/coherencefilter_version5b')
    
    % Initialize image data structure and pixelize (convert from nm to pix) setting values
    ims = initImgData(imfilei); % initialize the image structure, which generates a grayscale and some simple thresholded stuff to start with
    [settings, ims] = pix_settings(settings,ims);
    if settings.figSave
        mkdir(ims.imNamePath);  % make the directory to save the results figures
    end
    
    % Run the filter bank at the current settings
    auto_handles.ims = ims;
    auto_handles.settings = settings;
    auto_handles = main_filter(auto_handles);
    ims = auto_handles.ims;
    
    % Stitch fiber segments and calculate length
    ims = StitchFibers2(ims,settings);
    ims = FiberLengths(ims,settings);
    
    % Sample fiber widths along each fiber
    ims = FiberWidths(ims,settings);
    FiberData = [[ims.Fibers(:).Length]', [ims.Fibers(:).Width]', [ims.Fibers(:).Length]'./[ims.Fibers(:).Width]'];
    DataMat = [DataMat; FiberData];
    
%     save([ims.imNamePath, '_FiberData.mat'],'FiberData')

end

csvwrite([dirPath, 'fiber_data.csv'],DataMat);

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
