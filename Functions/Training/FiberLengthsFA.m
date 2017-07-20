function FLD = FiberLengthsFA(fa_path,noEdge)

if exist('noEdge')~=1
    noEdge=1;
end

load(fa_path)
xy = imageData.xy;
step_nm = imageData.step_nm;
sizeX = imageData.sizeX;
sizeY = imageData.sizeY;
numFibs = length(xy);
FLD = zeros(numFibs,1);

for i = 1:length(xy)
    
    xy_i=xy{i};
    length_i = (size(xy{i},2)-1)* step_nm;
    
    if noEdge
        
        if (not(any( xy_i(1,:) < 0.05 * sizeX )) && ...
                not(any( xy_i(1,:) > 0.95 * sizeX )) && ...
                not(any( xy_i(2,:) < 0.05 * sizeY )) && ...
                not(any( xy_i(2,:) > 0.95 * sizeY )) ...
                )
            
            FLD(i) = length_i;
        end
    else
        FLD(i) = length_i;
    end
end

FLD(FLD==0) = [];

end
