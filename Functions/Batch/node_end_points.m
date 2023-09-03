function xl = node_end_points(dirPath, settings, b)

% Compile images from the directory, run the full processing/analysis stack
% specified by the current GUI parameters, save results to csv and save
% visualizations to folder with image file's name
% b is the number of pixels around edge of image to exclude for endpoint
% analysis

hwaitdir = waitbar(0,'Running Directory...');

imdir = CompileImgs(dirPath);
numIms = length(imdir);
mkdir(fullfile(dirPath, "ep_results"));

% Initialize the Cell for the csv file
xl = cell(length(imdir)+1,3);
xl{1,1} = 'Image Name';
xl{1,2} = 'Endpoint Density (um^-2)';
% xl{1,3} = 'Node Density (um^-2)';

for i = 1:numIms
    
    waitbar(i/numIms,hwaitdir,['Processing ', imdir(i).name]);
    
    imfilei = imdir(i).path;
%     addpath(genpath(pwd));
    
    % Initialize image data structure
    ims = initImgData(imfilei);
    ims.settings = settings;
    ims = pix_settings(ims);
    
    % Run the filter bank at the current settings
    handles.ims=ims;
    handles = main_filter(handles);
    ims=handles.ims;
    
    % Generate branch and endpoints, calculate metrics
    [m,n] = size(ims.skelTrim);
    ims.branchpoints = bwmorph(ims.skelTrim, 'branchpoints');
    ims.bp_dens = sum(sum(ims.branchpoints)) / (m*ims.nmPix/1000*n*ims.nmPix/1000);
    endpoints_prelim = bwmorph(ims.skelTrim, 'endpoints');
    ims.endpoints = padarray(endpoints_prelim(1+b:m-b, 1+b:n-b), [b, b]); % remove spurious endpoints from image edges
    ims.ep_dens = sum(sum(ims.endpoints)) / (m*ims.nmPix/1000*n*ims.nmPix/1000);
    
    % Generate image showing branch and endpoints
    ims.bp_ep_mask = ims.img;
    % ims.bp_ep_mask = imoverlay(ims.bp_ep_mask, bwmorph(ims.branchpoints, 'dilate', 2), 'blue');
    ims.bp_ep_mask = imoverlay(ims.bp_ep_mask, bwmorph(ims.endpoints, 'dilate', 2), 'cyan');
    
    imwrite(ims.bp_ep_mask, fullfile(dirPath, "ep_results", ims.imName+"_bp_ep.png"))

    % Write data to csv cell
    xl{i+1,1} = imdir(i).name;
    xl{i+1,2} = ims.ep_dens;
    % xl{i+1,3} = ims.bp_dens;
    
end

savePath = fullfile(dirPath, "ep_results", "ep_results.csv");
cell2csv(savePath, xl, ',', 1999, '.');
close(hwaitdir)

end

