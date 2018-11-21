function ims = ODist_noVec(ims)

% Calculate the orientation distribution of pixels with valid angles
% (mask with skelTrim)

% Add fields to ims:

angleStep = 2;  % 2 degree discretization
stepNum = round(180/angleStep);
angleStep = 180/stepNum;

ang_list = ims.AngMap(ims.skelTrim)';  % These are from -90 to 90
rad_list = deg2rad(ang_list);

% Calculate for 0-pi range
n = histc(rad_list, linspace(-pi/2, pi/2, stepNum+1));
n(1) = n(1) + n(end);
n(end) = []; % remove orientation = pi values

% Reflect through the origin to the full 360 deg range
step = pi*angleStep/180; % angle step in rad
centers_polar = - step/2 + linspace(-pi/2, 3*pi/2, 2*stepNum+1);
centers_deg = 180*centers_polar/pi; % recalculate into deg in case of saving to a text file
n = [n(end), n, n];

ims.ODist.n = n;
ims.ODist.centers_polar = centers_polar;
ims.ODist.centers_deg = centers_deg;
ims.ODist.director = pi/2;
ims.ODist.S = 2*mean(cos(rad_list-pi/2).^2)-1;

end