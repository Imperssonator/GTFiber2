function [umFrames, sfull, smod, BETA] = op2d_am(ims, settings)

% Hard Coded
% -----------
fmin = 0;
fstep = settings.frameStep;
gridstep = settings.gridStep;
figSwitch = settings.figSwitch;
figSave = settings.figSave;
% -----------

angMap = ims.AngMap;
bw = ims.skelTrim;
nmPix = ims.nmPix;

[m,n] = size(angMap);
% disp(m); disp(n)
smallDim = min(m,n);
fmax = floor((smallDim-1)/2);

% Fix inputs to be round and non-negative
if fmin<0
    fmin = 0;
elseif fmin>fmax
    fmin = 0;
end
fmin = round(fmin);
fmax = round(fmax);
fstep = ceil(fstep);
gridstep = ceil(gridstep);
midFramei = ceil(m/2);
midFramej = ceil(n/2);

% Build frame size list, 'frames'
frames = (fmin:fstep:fmax);
if frames(end) ~= fmax
    frames = ceil([frames, fmax]);
end
numFrames = length(frames);

% for each frame size, sample angMap and get an average S2D
sfull = zeros(size(frames));

hwait2d = waitbar(0,'Calculating Order Parameter...');

for f = 1:numFrames
    fhw = frames(f);    % frame half width
    centi = (1+fhw:gridstep:m-fhw); % Coordinates of frame centers
    centj = (1+fhw:gridstep:n-fhw);
    midCenti = mean(centi);   % Shift the grid so that it is symmetrix about the center of the image
    midCentj = mean(centj);
    centShifti = midFramei-midCenti;
    centShiftj = midFramej-midCentj;
    centi = ceil(centi+centShifti);
    centj = ceil(centj+centShiftj);
    [fci, fcj] = meshgrid(centi,centj);
    [fm, fn] = size(fci);
    Sf = zeros(fm,fn);
    % Extract frames from angle map and calculate S2D for each
    for ii = 1:fm
        for jj = 1:fn
            curFrameAM = angMap(fci(ii,jj)-fhw:fci(ii,jj)+fhw,fcj(ii,jj)-fhw:fcj(ii,jj)+fhw);
            curFrameBW = bw(fci(ii,jj)-fhw:fci(ii,jj)+fhw,fcj(ii,jj)-fhw:fcj(ii,jj)+fhw);
            Sf(ii,jj) = angleFrameS2D(curFrameAM,curFrameBW);
        end
    end
    sfull(f) = nanmean(nanmean(Sf));    % Average S2D, excluding NaN values (from frames that had no fibers)
    waitbar(f/numFrames,hwait2d);
end

close(hwait2d)

% Convert frame sizes to microns and either make a figure or do the model
% fit to send back to the csv file if using runDir
umFrames = (frames.*2+1).*nmPix/1000;   % Frames in units of microns

% Remove NaN values from sfull if that happened by accident
umFrames = umFrames(~isnan(sfull));
sfull = sfull(~isnan(sfull));

if figSwitch || figSave
    [smod, BETA] = plotS2Dint(umFrames,sfull,ims,figSave);
else
%     BETA = nlinfit(umFrames,sfull,@(beta,x) beta(1)+(1-beta(1)).*exp(-x./beta(2)),[0.5, 1]);
    BETA = lsqnonlin(@(B) B(1)+(1-B(1)).*exp(-umFrames./B(2)) - sfull,[0.5,1],[0 1E-3],[1 Inf]);
    smod = BETA(1)+(1-BETA(1)).*exp(-umFrames./BETA(2));
end

end
    
function S = angleFrameS2D(curFrameAM,curFrameBW)

% 2D order parameter algorithm

orientList = curFrameAM(~~curFrameBW);
if isempty(orientList)
    S = NaN;
    return
end

olx = cosd(orientList);
oly = -sind(orientList);
vect = [olx'; oly'];
A = sum(vect(1,:).^2);      
B = sum(prod(vect));
N = size(vect, 2);
S = sqrt((2*A-N)^2+4*B^2)/N;

end

function [smod,BETA] = plotS2Dint(frames,sfull,ims,figSave)

% Use a finer frame size to plot the model (smoother)
fineFrames = linspace(frames(1),frames(end),500);

% BETA = nlinfit(frames,sfull,@(beta,x) beta(1)+(1-beta(1)).*exp(-x./beta(2)),[0.5, 1]);
BETA = lsqnonlin(@(B) B(1)+(1-B(1)).*exp(-frames./B(2)) - sfull,[0.5,1],[0 1E-3],[1 Inf]);
fitY = BETA(1)+(1-BETA(1)).*exp(-fineFrames./BETA(2));  % Use fineFrames to make figure look good
smod = BETA(1)+(1-BETA(1)).*exp(-frames./BETA(2));  % Use regular frames so that data lines up in csv

% Hard coded figure settings that look nice
if ispc
    font=24;
    flfont=22;
    flpos=[0.6, 0.855];
else
    font=30;
    flfont=28;
    flpos=[0.64, 0.87];
end
marker = 9;
markerline = 1.25;
line = 2.2;
lenscale = 1000;
edgedark = 0;
edgewidth = 0.75;

f=figure;
f.Position = [390   143   718   655];
ax = gca;
hold(ax,'on');
p1 = plot(ax,frames,sfull,'ok'); %, ...
p2 = plot(ax,fineFrames,fitY,'-b');
%     expt.dev(device).AFM(afm).xS/lenscale,expt.dev(device).AFM(afm).S_rand,'-r');
xlabel('Frame Size (µm)');
ylabel('{\itS}_{2D}')
ax.FontSize = font;
ax.Box = 'on';
ax.LineWidth = 0.75;
ax.PlotBoxAspectRatio = [1 1 1];
ax.YLim = [0 1];
ax.TickLength = [0.025, 0.025];
p1.LineWidth = markerline;
p1.MarkerSize = marker;
p2.LineWidth = line;
p2.Color = [0 0 1];
% p(2).LineWidth = line; p(2).Color = [1, 0, 0];

htex = text('Units', 'normalized', 'Position', flpos, ...
    'BackgroundColor', [1 1 1], ...
    'String', {['{\it\lambda}_{C}= ', num2str(BETA(2)*1000,4), ' nm'],...
               ['{\itS}_{full}= ', num2str(BETA(1),2)]},...
    'FontSize', flfont,...
    'EdgeColor', edgedark*[1 1 1],...
    'LineWidth', edgewidth);

if figSave
    hgexport(f, [ims.figSavePath, '_OP2D', '.tif'],  ...
        hgexport('factorystyle'), 'Format', 'tiff');
    close(f)
end

end
