function handles = main_filter(handles)

settings = handles.ims.settings;
ims = handles.ims;

% Run coherence filter
hwait = waitbar(0,'Diffusion Filter...');
Options = settings.Options;
switch settings.invert
    case 0
        [ims.CEDgray, ims.v1x, ims.v1y] ...
            = CoherenceFilter(imadjust(ims.gray),Options);
        ims.CEDgray = mat2gray(ims.CEDgray);
    case 1
        [ims.CEDgray, ims.v1x, ims.v1y] ...
            = CoherenceFilter(imcomplement(ims.gray),Options);
        ims.CEDgray = mat2gray(ims.CEDgray);
end
handles=imshowGT(ims.CEDgray,handles,'img_axes');


% Run Top Hat Filter
waitbar(0.5,hwait,'Top Hat Filter...');
ims.CEDtophat = imadjust(imtophat(ims.CEDgray,strel('disk',settings.thpix)));
handles=imshowGT(ims.CEDtophat,handles,'img_axes');


% Threshold and Clean
waitbar(0.7,hwait,'Threshold and Clean...');
switch settings.threshMethod
    case 1
        ims.CEDbw = imbinarize(ims.CEDtophat,'adaptive');
    case 2
        ims.CEDbw = im2bw(ims.CEDtophat,settings.globalThresh);
end

handles=imshowGT(ims.CEDbw,handles,'img_axes');
ims.CEDclean = bwareaopen(ims.CEDbw,settings.noisepix);
% ims.cleanRP = regionprops(ims.CEDclean,'Solidity','Eccentricity');
% not_particles = ~([ims.cleanRP(:).Eccentricity]<0.95 & [ims.cleanRP(:).Solidity]>0.8);
% temp_label = bwlabel(ims.CEDclean);
% ims.CEDclean = MultiEquiv(temp_label,find(not_particles));
% ims.CEDclean = imclose(ims.CEDclean,strel('disk',1));
handles=imshowGT(ims.CEDclean,handles,'img_axes');


% Skeletonize
waitbar(0.8,hwait,'Skeletonization...');
ims.skel = bwmorph(ims.CEDclean,'thin',Inf);
handles=imshowGT(ims.skel,handles,'img_axes');
ims = CleanSkeleton(ims);
handles=imshowGT(ims.skelTrim,handles,'img_axes');


% Generate Angle Map by getting new angles from CED
waitbar(0.9,hwait,'Recovering Orientations...');
ims.AngMap = atand(ims.v1x./-ims.v1y);

handles.ims = ims;
% save('filter_debug','ims')
close(hwait)

end
