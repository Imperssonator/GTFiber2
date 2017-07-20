function [umFrames, sfull, smod, BETA] = gif_op2d_am(ims, settings)

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

% Initialize gif file
gif_file = [ims.imNamePath, '_op2d.gif'];
gif_fiber_imgs = cell(numFrames+1,1);      % Where to store the fiber frame images

for f = 1:numFrames
    fhw = frames(f);    % frame half width
    
    centi = (1+fhw:gridstep:m-fhw); % Coordinates of frame centers, from frame half width to end - frame half width, by gridstep
    centj = (1+fhw:gridstep:n-fhw);
    midCenti = mean(centi);   % Shift the grid so that it is symmetrix about the center of the image
    midCentj = mean(centj);
    centShifti = midFramei-midCenti;  % difference between the actual middle of the image and the middle-most frame center
    centShiftj = midFramej-midCentj;
    centi = ceil(centi+centShifti);     % shift all the frame centers
    centj = ceil(centj+centShiftj);
    [fci, fcj] = meshgrid(centi,centj); % Make the meshgrid of the frame centers
    [fm, fn] = size(fci);
    
    gif_fiber_imgs{f} = bw(midFramei-fhw:midFramei+fhw,midFramej-fhw:midFramej+fhw);

    Sf = zeros(fm,fn);  % Where we store all of the S2D's extracted from each frame at this frame size
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

% Convert black and whites into gray
gif_fiber_imgs{numFrames+1} = bw;
gif_fiber_imgs = cellfun(@(x) uint8(x(:,:,[1 1 1]).*254+1),gif_fiber_imgs,'UniformOutput',false); % convert to "RGB" images

close(hwait2d)

% Convert frame sizes to microns and either make a figure or do the model
% fit to send back to the csv file if using runDir
umFrames = (frames.*2+1).*nmPix/1000;
if figSwitch || figSave
    [smod, BETA,gif_plot_imgs] = plotS2Dint(umFrames,sfull,ims,settings);
else
    BETA = nlinfit(umFrames,sfull,@(beta,x) beta(1)+(1-beta(1)).*exp(-x./beta(2)),[0.5, 1]);
    smod = BETA(1)+(1-BETA(1)).*exp(-umFrames./BETA(2));
end

% Now we have gif_fiber_imgs and gif_plot_imgs, and we simply concatenate
% them along the 2nd dimension and turn that sequence of image into a gif

plot_height = size(gif_plot_imgs{1},1);
gif_fiber_imgs = cellfun(@(IM) imresize(IM,[plot_height,plot_height]),gif_fiber_imgs,'UniformOutput',false);

% fiber_sizes = cellfun(@(x) size(x),gif_fiber_imgs,'UniformOutput',false);
% plot_sizes = cellfun(@(x) size(x),gif_plot_imgs,'UniformOutput',false);
% save('gif_debug')
% disp('saved')

gif_cat_imgs = cellfun(@(fib,plot) cat(2,fib,plot),gif_fiber_imgs,gif_plot_imgs,'UniformOutput',false);
% gif_cat_inds = cellfun(@(IM) rgb2ind(IM,256),gif_cat_imgs,'UniformOutput',false);
% disp('concat')

[gif_frame_1,cmap1] = rgb2ind(gif_cat_imgs{1},256);
imwrite(gif_frame_1,cmap1,gif_file,'gif','LoopCount',1,'DelayTime',settings.plotDelay);

for i = 2:numFrames
    [gif_frame_i,cmapi] = rgb2ind(gif_cat_imgs{i},256);
    imwrite(gif_frame_i,cmapi,gif_file,'gif','WriteMode','append','DelayTime',settings.plotDelay);
end

[gif_frame_end,cmap_end] = rgb2ind(gif_cat_imgs{end},256);
imwrite(gif_frame_end,cmap_end,gif_file,'gif','WriteMode','append','DelayTime',settings.plotFinal);
% disp('image saved')

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

function [smod,BETA,gif_plot_imgs] = plotS2Dint(frames,sfull,ims,settings)

numFrames = length(frames);
gif_plot_imgs = cell(numFrames+1,1);      % Where to store the fiber frame images


% Use a finer frame size to plot the model (smoother)
fineFrames = linspace(frames(1),frames(end),500);

BETA = nlinfit(frames,sfull,@(beta,x) beta(1)+(1-beta(1)).*exp(-x./beta(2)),[0.5, 1]);
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

% Start making figure

for i = 1:numFrames+1
    f_i=figure;
    f_i.Position = [390   143   718   655];
    ax = gca;
    hold(ax,'on');
    
    j=i;    % This is hacky but it's so we only plot 1:numFrames when we're making the final plot
    if i == numFrames+1
        p2 = plot(ax,fineFrames,fitY,'-b');
        j = i-1;
    end
    p1 = plot(ax,frames(1:j),sfull(1:j),'ok');  % Plot from 1 to the current frame
    
    
    xlabel('Frame Size (µm)');
    ylabel('{\itS}_{2D}')
    ax.FontSize = font;
    ax.Box = 'on';
    ax.LineWidth = 0.75;
    ax.PlotBoxAspectRatio = [1 1 1];
    ax.XLim = [0 ceil(settings.nmWid/1000)];
    ax.YLim = [0 1];
    ax.TickLength = [0.025, 0.025];
    p1.LineWidth = markerline;
    p1.MarkerSize = marker;
    
    if i == numFrames+1
        p2.LineWidth = line;
        p2.Color = [0 0 1];
    end
    % p(2).LineWidth = line; p(2).Color = [1, 0, 0];

    htex = text('Units', 'normalized', 'Position', flpos, ...
        'BackgroundColor', [1 1 1], ...
        'String', {['{\it\lambda}_{C}= ', num2str(BETA(2)*1000,4), ' nm'],...
                   ['{\itS}_{full}= ', num2str(BETA(1),2)]},...
        'FontSize', flfont,...
        'EdgeColor', edgedark*[1 1 1],...
        'LineWidth', edgewidth);

    gif_plot_imgs{i} = frame2im(getframe(f_i));     % Capture the plot and store it
    close(f_i)
end

end
