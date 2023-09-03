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
    CurIms(p).path = fullfile(FolderPath, CurIms(p).name);   % prepend the folder path to the image names
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
