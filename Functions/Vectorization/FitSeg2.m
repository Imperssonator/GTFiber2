function [xy_cur,curv,MSE,L] = FitSeg2(ims,segNum,params)

% This function is meant to stretch a contour until it hits the edges of an
% image

% Params:
% alpha
% beta
% gamma
% k1
% k2
% fiberIntensity
% Iterations

% meanInt = mean(ims.gray(:));

xy_init=ims.fibSegs(segNum).xy;

% fibSegImg = uint8(ims.CEDgray.*255);
% fiberStep = ims.settings.fiberStep;
% bgInt = mean(fibSegImg(~ims.CEDclean));
% fiberIntensity = mean(fibSegImg(ims.SegLabels==segNum));

fibSegImg = uint8(ims.CEDgray.*255);
fiberStep = ims.settings.fiberStep;
bgInt = mean(fibSegImg(~ims.CEDclean));
fiberIntensity = mean(fibSegImg(ims.SegLabels==segNum));


[gradX, gradY] = gradient2Dx2(int16(fibSegImg));

% Short names for fiber tracking parameters and data
xy_cur = distributePoints(xy_init,fiberStep);

a = params(1); %0;
b = params(2); %5;
g = params(3); %10;
k1 = params(4); %20;
k2 = params(5); %10;
% fiberIntensity = params(6); %255;

MSE = zeros(params(7),1);
L = zeros(params(7),1);

for k = 1:params(7)
    disp(k)
    % Construct a matrix M for the internal energy contribution
    n = length(xy_cur);
    
    m1 = [3   2 1:n-2         ;
        2     1:n-1         ;
        1:n           ;
        2:n   n-1     ;
        3:n   n-1 n-2];
    
    m2 = cumsum(ones(5,n), 2);
    
    col = [b; -a-4*b; 2*a+6*b+g; -a-4*b; b];
    m3 = col(:, ones(n,1));
    
    m3(:,1:2) = [b       0        ;
        -2*b    -a-2*b   ;
        a+2*b+g 2*a+5*b+g;
        -a-2*b  -a-4*b   ;
        b       b       ];
    
    m3(:,end-1:end) = [b         b      ;
        -a-4*b    -a-2*b ;
        2*a+5*b+g a+2*b+g;
        -a-2*b    -2*b   ;
        0         b     ];
    
    M = sparse(m1, m2, m3, n, n);
    
    % Vector field of external energy
    % (interp2D MATLAB function requests double input, which consume too
    % much memory for big AFM images)
    % gradX and gradY are doubled! (see gradient2Dx2)
    vfx = interp2D(gradX, xy_cur(1,:), xy_cur(2,:), 'cubic')/2;
    vfy = interp2D(gradY, xy_cur(1,:), xy_cur(2,:), 'cubic')/2;
    vf = k1*[vfx; vfy]/(fiberIntensity-bgInt);
    
    % Normalized vectors of fiber ends
    v_start = xy_cur(:,1) - xy_cur(:,2);
    v_start = v_start/sqrt(sum(v_start.^2));
    v_end = xy_cur(:,end) - xy_cur(:,end-1);
    v_end = v_end/sqrt(sum(v_end.^2));
    
    % Image intensities at the ends
    int_start = interp2D(fibSegImg, xy_cur(1,1), xy_cur(2,1), 'cubic');
    int_end = interp2D(fibSegImg, xy_cur(1,end), xy_cur(2,end), 'cubic');
    
    disp('----')
    extend_start = k2*v_start*(int_start-bgInt)/(fiberIntensity-bgInt)
    extend_end = k2*v_end*(int_end-bgInt)/(fiberIntensity-bgInt)
    disp('----')
    
    % Full external energy
    vf(:,1) = vf(:,1) + extend_start;
    vf(:,end) = vf(:,end) + extend_end;
    
    % Next position
    xy_next = distributePoints((g*xy_cur+vf)/M, fiberStep);
    if length(xy_next)==length(xy_cur)
        MSE(k)=mean(sum((xy_next-xy_cur).^2,1));
    else
        MSE(k)=1;
    end
    L(k) = length(xy_next);
    
%     Stop if equilibrated
%     if MSE(k)<0.005 %k>10 && ((MSE(k)-MSE(k-1))^2/MSE(k))<0.00001
%         %         disp(k)
%         break
%     end
    
    xy_cur = xy_next;
    %     % Stop if at edges
    %     if at_edge(xy_cur(:,1),size(fibSegImg)) && at_edge(xy_cur(:,end),size(fibSegImg))
    %         break
    %     end
end

if size(xy_cur,2)>=3
    dxy=(xy_cur(:,3:end)-xy_cur(:,1:end-2))./(2*fiberStep);
    ddxy=(xy_cur(:,3:end)-2.*xy_cur(:,2:end-1)+xy_cur(:,1:end-2))./(fiberStep)^2;
    curv=(dxy(1,:).*ddxy(2,:)-dxy(2,:).*ddxy(1,:))./sum(dxy.^2,1).^(1.5);
    curv=[0,abs(curv),0];
else
    curv = zeros(1,size(xy_cur,2));
end

MSE(MSE==0)=[];
L(L==0)=[];
% disp(MSE(k))

end

function out = at_edge(xy,imSize)

if xy(1) <= 0.01*imSize(2) || xy(1) >= 0.99*imSize(2) ...
        || xy(2) <= 0.005*imSize(1) || xy(2) >= 0.995*imSize(1)
    out=1;
else
    out=0;
end


end


