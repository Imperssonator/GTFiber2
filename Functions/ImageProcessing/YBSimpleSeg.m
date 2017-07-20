function bw = YBSimpleSeg(gray)

%YBSeg Yanowitz-Bruckstein image segmentation with fiber unentanglement
%   SP is the structure path
%   File is the file path from current active dir
%   Dim is the image dimension in nm

edges = edge(gray,'canny');         % Apply a Canny edge finder - 1's at the edges, 0's elsewhere
grayDouble = double(gray);                     % Turn the grey image into double prec.
edgesDouble = double(edges);                     % Turn the edge image into double prec.
initThresh = grayDouble.*edgesDouble;                        % Fill in the grey values of the edge pixels in a new image file                      

threshSurf = YBiter(initThresh);                    % Perform Yanowitz-Bruckstein surface interpolation to create threshold surface from edge gray values

bw = gray>threshSurf;                          % Segment the image; pixels above threshold surface are white, if not, black

end

function Vf = YBiter(V0)

%YBiter Yanowitz/Bruckstein surface interpolation

w = 1;
[m,n] = size(V0);
InitVal = mean(V0(V0~=0)); % Start non-edges as the mean of the edges because 0's aint' working
Vupdate = V0==0;           % A logical array of pixels to update on each iteration
Vupdate(1,:) = 0;
Vupdate(:,1) = 0;
Vupdate(m,:) = 0;
Vupdate(:,n) = 0;

V0(V0==0)=InitVal;         % Put the average edge values in the non-edge cells

Vnew = zeros(m,n);         % Initialize the updated threshold surface
Vold = V0;

iter = 0;
maxiter = 40;

% tic
while iter<maxiter
    iter = iter+1;
%     disp(iter)
    
    Lap = del2(Vold);
    Vnew = Vold + Vupdate .* (w .* Lap);
    Vold = Vnew;
    
end
% toc
Vf = Vnew;

end