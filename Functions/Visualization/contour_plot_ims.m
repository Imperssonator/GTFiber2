function [him, hs, hax, hf, hc] = contour_plot_ims(ims,segorfib,segNums)

switch segorfib
    case 'seg'
        
        if exist('segNums')~=1
            segNums = (1:length(ims.fibSegs));
        end
        
        xy_sel = {ims.fibSegs(segNums).xy};
        
        [him, hp, ho, hax, hf] = ...
            contour_plot(...
            ims.gray,...
            xy_sel);
        
    case 'fib'
        
        if exist('segNums')~=1
            segNums = (1:length(ims.Fibers));
        end
        
        xy_sel = {ims.Fibers(segNums).xy};
        
        [him, hp, ho, hax, hf] = ...
            contour_plot(...
            ims.gray,...
            xy_sel);
end

end