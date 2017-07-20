function ims = RemoveBadWidths(ims)

% Start with Initial segs
segsWidthSpec = ims.segsInit;

for i = 1:length(ims.fibSegs)
    
    if ims.fibSegs(i).width < ims.settings.minWidthnm || ims.fibSegs(i).width > ims.settings.maxWidthnm
        segsWidthSpec(ims.fibSegs(i).sortPixInds) = 0;
    end

end

ims.segsInit = segsWidthSpec;

% And then we restart the whole process!

end