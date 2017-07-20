load('AC_example')

% Initialize 3 points along fiber
xy_init=[300 325 350;...
         37 32 32];
fiberStep=3;
xy_init=xy_init;

% Distribute 
xy = distributePoints(xy_init,fiberStep);
xy_dist=xy;

f=figure;ax=gca;
hold on
h_im=imshow(ims.gray);
ax.XLim=[290 365];
ax.YLim=[21 44];
h_dist=plot(ax,xy(1,:),xy(2,:),'ob');
% h_init=plot(ax,xy_init(1,:),xy_init(2,:),'or');
% h_init.MarkerSize=15;
% h_init.LineWidth=2;

h_im.AlphaData=0.4;
h_dist.MarkerSize=10;
h_dist.LineWidth=2;

% Run one iteration
[gradX, gradY] = gradient2Dx2(int16(ims.gray));
hq = quiver(ax,gradX,gradY);
hq.AutoScaleFactor=1.5;
hq.LineWidth=1;

a = 0;
b = 5;
g = 10;
k1 = 20;
k2 = 10;

fiberIntensity = 255;
fibImg=ims.gray;

for k = 1:20
    
    n = length(xy);
    
    m1 = [3   2 1:n-2         ;
        2     1:n-1         ;
        1:n           ;
        2:n   n-1     ;
        3:n   n-1 n-2];
    
    m2 = cumsum(ones(5,n), 2);
    
    col = [b; -a-4*b; 2*a+6*b+g; -a-4*b; b];
    
    m3 = col(:, ones(n,1));
    m3(:,1:2) = ...
    [b       0        ;
     -2*b    -a-2*b   ;
     a+2*b+g 2*a+5*b+g;
     -a-2*b  -a-4*b   ;
     b       b       ];
    
    m3(:,end-1:end) = ...
       [b         b      ;
        -a-4*b    -a-2*b ;
        2*a+5*b+g a+2*b+g;
        -a-2*b    -2*b   ;
        0         b     ];
    
    M = sparse(m1, m2, m3, n, n);
    
    % Vector field of external energy
    % (interp2 MATLAB function requests double input, which consume too
    % much memory for big AFM images)
    % gradX and gradY are doubled! (see gradient2Dx2)
    vfx = interp2D(gradX, xy(1,:), xy(2,:), 'cubic')/2;
    vfy = interp2D(gradY, xy(1,:), xy(2,:), 'cubic')/2;
    vf = k1*[vfx; vfy]/fiberIntensity;
    
    % Normalized vectors of fiber ends
    v_start = xy(:,1) - xy(:,2);
    v_start = v_start/sqrt(sum(v_start.^2));
    v_end = xy(:,end) - xy(:,end-1);
    v_end = v_end/sqrt(sum(v_end.^2));
    
    % Image intensities at the ends
    int_start = interp2D(fibImg, xy(1,1), xy(2,1), 'cubic');
    int_end = interp2D(fibImg, xy(1,end), xy(2,end), 'cubic');
    
    % Full external energy
    vf(:,1) = vf(:,1) + k2*v_start*int_start/fiberIntensity;
    vf(:,end) = vf(:,end) + k2*v_end*int_end/fiberIntensity;
    
    xy = distributePoints((g*xy+vf)/M, fiberStep);
    h_new=plot(ax,xy(1,:),xy(2,:),'.m');
    h_new.MarkerSize=10;h_new.LineWidth=2;
    
end

h_new=plot(ax,xy(1,:),xy(2,:),'om');
h_new.MarkerSize=10;h_new.LineWidth=2;

ax.Position=[0 0 1 1];
f.Position=[353 449 953 293];