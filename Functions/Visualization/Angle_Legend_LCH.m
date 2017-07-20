function legfig = Angle_Legend_LCH(pos)

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

% Bow the color distribution a bit because cyan is too powerful
pow = 0.5;
% Move the angles forward so cyan is at the end
rot = 35;
Angs = atan2d(Y,X)*2;
Angs = Angs-rot;
Angs(Angs<0) = Angs(Angs<0)+360;
H = (Angs).^pow./(360^(pow-1));
C = ones(size(X))*100;
L = ones(size(X));
A = ones(size(X));

% Make anything outside the radius of the color wheel transparent
Clear = R>rr;
Color = R<=rr;
L(Clear) = 0;
L(Color) = 55;
A(Clear) = 0;
A(Color) = 1;

% Convert to rgb
LCH = cat(3,L,C,H);
rgb = colorspace('LCH->RGB',LCH);

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