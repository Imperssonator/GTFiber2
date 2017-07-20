function ims = calc_orient_dist(ims)

% Calculate the orientation distribution of segment vectors from all of the
% fitted fibers in ims

% Add fields to ims:

angleStep = 5;  % 5 degree discretization
stepNum = round(180/angleStep);
angleStep = 180/stepNum;

xy = {ims.Fibers(:).xy};  % This is how FiberApp handles the xy data
vect = cellfun(@get_vectors, xy, 'UniformOutput', false);
vect = [vect{:}]; % unite data of all fibrils

% Calculate 2D order parameter (S)
A = sum(vect(1,:).^2);
B = sum(prod(vect));
N = size(vect, 2);
S = sqrt((2*A-N)^2+4*B^2)/N;

% Calculate orientation distribution
vect(:, vect(2,:)>0) = - vect(:, vect(2,:)>0); % Turn the coordinate system from informatics into geometrical 
orientation = acos(vect(1,:));

% Calculate Average Orientation
J = zeros(length(orientation),1,4);
J(:,1,1) = cos(orientation).^2;
J(:,1,2) = cos(orientation).*sin(orientation);
J(:,1,3) = cos(orientation).*sin(orientation);
J(:,1,4) = sin(orientation).^2;
MeanJ = mean(J,1);
MeanOrient = RecoverAnglesOD(MeanJ);
if MeanOrient<0
    MeanOrient = MeanOrient + 180;
end

% Calculate for 0-pi range
n = histc(orientation, linspace(0, pi, stepNum+1));
n(1) = n(1) + n(end);
n(end) = []; % remove orientation = pi values

% Reflect through the origin to the full 360 deg range
step = pi*angleStep/180; % angle step in rad
centers_polar = - step/2 + linspace(0, 2*pi, 2*stepNum+1);
centers_deg = 180*centers_polar/pi; % recalculate into deg in case of saving to a text file
n = [n(end), n, n];

ims.ODist.n = n;
ims.ODist.centers_polar = centers_polar;
ims.ODist.centers_deg = centers_deg;
ims.ODist.director = MeanOrient;

end

function v = get_vectors(xy)
v = diff(xy, 1, 2);
l = sqrt(sum(v.^2));
v = v./[l; l];

end

function [AngMap, Dirs] = RecoverAnglesOD(J)

% disp('Finding Angles...')

[m, n] = size(J(:,:,1));
JE = zeros(m,n,2);
JV = zeros(m,n,4);
for i = 1:m
    for j = 1:n
        [V,D] = eig(reshape(J(i,j,:),[2 2]));
        JE(i,j,:) = diag(D);
        JV(i,j,:) = V(:);
    end
end
JESort = sort(JE,3,'descend');
JVSort = zeros(size(JV));
for i = 1:m
    for j = 1:n
        if JESort(i,j,1)~=JE(i,j,1)
            JVSort(i,j,:) = [JV(i,j,3) JV(i,j,4) JV(i,j,1) JV(i,j,2)];
        else
            JVSort(i,j,:) = JV(i,j,:);
        end
    end
end

% Coher = ((JESort(:,:,1)-JESort(:,:,2))./(JESort(:,:,1)+JESort(:,:,2))).^2;
Dirs = JVSort(:,:,1:2);
AngMap = atan2d(JVSort(:,:,2),JVSort(:,:,1));      % last evec is 'coherence orientation'


end

