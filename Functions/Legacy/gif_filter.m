function gif_filter(ims,settings)

Options = settings.Options;
gif_file = [ims.imNamePath, '.gif'];

% Write initial gray image
[curFrame,c_map] = gray2ind(ims.gray,256);
imwrite(curFrame,c_map,gif_file,'gif','LoopCount',3,'DelayTime',settings.initDelay);

% Run coherence filter
hwait = waitbar(0,'Diffusion Filter...');
[ims.CEDgray, ims.v1x, ims.v1y, CED_gif] ...
    = CoherenceFilter(ims.gray,Options);
ims.CEDgray = mat2gray(ims.CEDgray);

% Write steps of CED
for i = 1:length(CED_gif)-1
    [curFrame, c_map] = gray2ind(CED_gif{i},256);
    imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDStepDelay);
end
i = length(CED_gif);
[curFrame, c_map] = gray2ind(CED_gif{i},256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);

% Run Top Hat Filter
waitbar(0.5,hwait,'Top Hat Filter...');
ims.CEDtophat_init = imtophat(ims.CEDgray,strel('disk',settings.thpix));
ims.CEDtophat = imadjust(ims.CEDtophat_init);

% Put Top Hat Filter in the Gif
[curFrame, c_map] = gray2ind(ims.CEDtophat_init,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);
[curFrame, c_map] = gray2ind(ims.CEDtophat,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);


% Threshold and Clean
waitbar(0.7,hwait,'Threshold and Clean...');
switch settings.threshMethod
    case 1
        ims.CEDbw = YBSimpleSeg(ims.CEDtophat);
        disp('used adaptive')
    case 2
        ims.CEDbw = im2bw(ims.CEDtophat,settings.globalThresh);
end

% Put Threshold in the Gif
[curFrame, c_map] = gray2ind(ims.CEDbw,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);

ims.CEDclean = bwareaopen(ims.CEDbw,settings.noisepix);

% Put Clean in the Gif
[curFrame, c_map] = gray2ind(ims.CEDclean,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);

% Skeletonize
waitbar(0.8,hwait,'Skeletonization...');
ims.skel = bwmorph(ims.CEDclean,'skel',Inf);

% Put Skeleton in the Gif
[curFrame, c_map] = gray2ind(ims.skel,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.CEDFinalDelay);

ims.skelTrim = cleanSkel(ims.skel,settings.maxBranchSize);
% ims.skelTrim = bwareaopen(ims.skelTrim,settings.maxStubLen);

% Put Trim Skeleton in the Gif
[curFrame, c_map] = gray2ind(ims.skelTrim,256);
imwrite(curFrame,c_map,gif_file,'gif','WriteMode','append','DelayTime',settings.skelDelay);

close(hwait)

end
