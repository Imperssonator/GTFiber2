function ims = pix_settings(ims)

% Conversions between nm and pixels
ims.nmWid = ims.settings.nmWid;
ims.pixWid = size(ims.img,2);
ims.nmPix = ims.nmWid/ims.pixWid;

% Get pixel values for filter options
ims.settings.Options.sigma = ims.settings.gaussnm/ims.nmPix;
ims.settings.Options.rho = ims.settings.rhonm/ims.nmPix;

% Top Hat
ims.settings.thpix = ceil(ims.settings.thnm/ims.nmPix);

% Noise Removal
ims.settings.noisepix = ceil(ims.settings.noisenm/ims.nmPix^2);

% Skeletonization
ims.settings.maxBranchSize = ceil(ims.settings.maxBranchSizenm/ims.nmPix);

% Match Search settings
ims.settings.fiberStep = ims.settings.fiberStep_nm / ims.nmPix;
ims.settings.searchLat = 2 * ims.settings.stitchGap;
ims.settings.searchLong = 2 * ims.settings.stitchGap;


end