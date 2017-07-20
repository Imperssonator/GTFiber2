function legfig = Angle_Legend(pos)

% Create a hemispherical color wheel as a legend for the colored angles in
% the angle color map

if exist('pos','var')~=1
    pos=[0 0];
end

% Set radius of hemisphere in pixels
rr = 220;
h = rr;
wl = -rr; wr = rr;

% Create a meshgrid, convert it to polar coordinates, and define HSV color
% values based on those coordinates.
[X, Y] = meshgrid((wl:wr),(h:-1:0));
R = sqrt((X.^2+Y.^2));
H = atan2d(Y,X)*2;
S = ones(size(X));
L = ones(size(X));
A = ones(size(X));

% Make anything outside the radius of the color wheel transparent
Clear = R>rr;
Color = R<=rr;
L(Color) = 0.6;
L(Clear) = 0;
A(Clear) = 0;
A(Color) = 1;

% Convert to rgb
HSL = cat(3,H,S,L);
rgb = colorspace('HSL->RGB',HSL);

% Add angle labels in degrees
FontSize=80;

rgbt1 = insertText(rgb,[-10 110],['180', char(176)],'BoxColor','black','TextColor','white','FontSize',FontSize,'BoxOpacity',0);
rgbt2 = insertText(rgbt1,[315 110],['0', char(176)],'BoxColor','black','TextColor','white','FontSize',FontSize,'BoxOpacity',0);
rgbt3 = insertText(rgbt2,[135 -20],['90', char(176)],'BoxColor','black','TextColor','white','FontSize',FontSize,'BoxOpacity',0);

Legend = cat(3,rgbt3,A);

legfig = figure;
him = imshow(rgbt3);
set(him,'AlphaData',A);
legfig.Position(1:2) = pos;
legfig.Children.Position = [0 0 1 1];
legfig.Position(3:4) = [2*rr, rr];

end