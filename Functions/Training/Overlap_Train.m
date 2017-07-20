if ispc
    im_path = 'C:\Users\npersson3\Google Drive\GTFiber Testing\SonAge Sol3 blade 2 V3.002 rot.tif';
    fa_path = 'C:\Users\npersson3\Google Drive\GTFiber Testing\SonAge Sol3 blade 2 V3.002 rot.mat';
else
    im_path = '/Users/Imperssonator/Google Drive/Kaylie/Fiber Growth master/For FiberApp/S2U0T20B2D0126.001.tif';
    fa_path = '/Users/Imperssonator/Google Drive/Kaylie/Fiber Growth master/For FiberApp/S2U0T20B2D0126.001a.mat';
end
nmWid=5000;
optimize_train_pt(im_path,fa_path,nmWid,'SonAgeOptim');


im_path = 'Example Images\Fig 5C 5000nm.tif';
fa_path = 'Example Images\Fig 5C 5000nm.mat';

optimize_train_pt(im_path,fa_path,nmWid,'MFUOptim')