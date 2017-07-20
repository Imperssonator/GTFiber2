function ims = RemoveShortFibers(ims)

FibersNew = ims.Fibers(1);

count = 0;
for i = 1:length(ims.Fibers)
    if ims.Fibers(i).length > ims.settings.minFibLen
        count = count+1;
        FibersNew(count) = ims.Fibers(i);
    end
end

rmfield(ims,'Fibers')
ims.Fibers = FibersNew;

% Build the FiberLabels matrix
FiberLabels = zeros(size(ims.SegLabels));
for i = 1:length(ims.Fibers)
    FiberLabels = FiberLabels ...
        + i .* MultiEquiv(ims.SegLabels,ims.Fibers(i).Fiber);
end

ims.FiberLabels = FiberLabels;

end