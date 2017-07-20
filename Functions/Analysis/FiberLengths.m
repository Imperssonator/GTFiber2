function ims = FiberLengths(ims,noEdge)

if exist('noEdge')~=1
    noEdge=1;
end

numFibs = length(ims.Fibers);
FLD = zeros(numFibs,1);

for i = 1:length(ims.Fibers)
    
    length_i = ims.Fibers(i).length;
    
    if noEdge
        if (not(any( ims.Fibers(i).xy(1,:) < 0.05 * size(ims.gray,2) )) && ...
                not(any( ims.Fibers(i).xy(1,:) > 0.95 * size(ims.gray,2) )) && ...
                not(any( ims.Fibers(i).xy(2,:) < 0.05 * size(ims.gray,1) )) && ...
                not(any( ims.Fibers(i).xy(2,:) > 0.95 * size(ims.gray,1) )) ...
                )
            
            FLD(i) = length_i;
        end
    else
        FLD(i) = length_i;
    end
end

FLD(FLD==0) = [];
ims.FLD = FLD;

end
