function xl = run_directory_CS(handles,dirPath,savePath,Mn_vec,PDI_vec,nruns)

%run_directory_CS
% Compile images from the directory, run the full processing/analysis stack
% specified by the current GUI parameters, then run the chain stacking
% simulation for each molecular weight and PDI specified, with nruns
% realizations of the chain stacking simulation for each
% 
% Mn_vec is nx1 list of molecular weights in Daltons (ex. 10000)
% PDI_vec is mx1 list of polydispersities (unitless) (ex. 1.5)
% 
% Total number of simulations run for each image is nruns * n * m, so total
% number of simulations run overall is numIms * nruns * n * m


hwaitdir = waitbar(0,'Running Directory...');

% Compile image files in folder
imdir = CompileImgs(dirPath);
numIms = length(imdir);
numMn = size(Mn_vec,1);
numPDI = size(PDI_vec,1);

% Initialize the Cell for the csv file
xl = cell(numIms*numMn*numPDI+1,14);

xl{1,1} = 'Image Name';
xl{1,2} = 'Sfull fit';
xl{1,3} = 'Correlation Length (nm)';
xl{1,4} = 'Average Orientation (degrees)';
xl{1,5} = 'Fiber Length Density (1/um)';
xl{1,6} = 'Mean Fiber Length (nm)';
xl{1,7} = 'Mean Fiber Width (nm)';
xl{1,8} = 'Mn (kDa)';
xl{1,9} = 'PDI';
xl{1,10} = 'Tie Chain Density (1/nm)';
xl{1,11} = 'TCD std. dev.';
xl{1,12} = 'Number of simulations';
xl{1,13} = 'Number of Fibers';
xl{1,14} = 'Number of Chains';

% Start count
count = 1;

for i = 1:numIms
    
    waitbar(i/numIms,hwaitdir,['Processing ', imdir(i).name]);
    
    imfilei = imdir(i).path;
    addpath(genpath(pwd));
    
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
    
    % Save figures if specified
    if ims.settings.figSave
        FiberVec_ACM(ims,1);
        ODist_plot(ims,1);
        plotS2D(ims,1);
    end
    
    for m = 1:numMn
        for p = 1:numPDI
            
            count = count+1;
            
            tcd_vec = [];
            for r = 1:nruns
                ims_temp = populate_chains(ims,Mn_vec(m),PDI_vec(p),0);
                ims_temp = tie_chain_density(ims_temp);
                tcd_vec = [tcd_vec; mean(ims_temp.tc_dist)];
            end
            
            % Write data to csv cell
            xl{count,1} = imdir(i).name;
            xl{count,2} = ims.op2d.Sfull;
            xl{count,3} = ims.op2d.decayLen;
            xl{count,4} = ims.ODist.director;
            xl{count,5} = ims.fibLengthDensity;
            xl{count,6} = mean(ims.FLD);
            xl{count,7} = mean(ims.FWD);
            xl{count,8} = Mn_vec(m);
            xl{count,9} = PDI_vec(p);
            xl{count,10} = mean(tcd_vec);
            xl{count,11} = std(tcd_vec);
            xl{count,12} = nruns;
            xl{count,13} = length(ims_temp.Fibers);
            xl{count,14} = size(ims_temp.Chains,1);
            
        end
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
