function [him, hs, hax, hf, hc] = curv_plot_ims(ims,segorfib,segNums)

switch segorfib
    case 'seg'
        
        if exist('segNums')~=1
            segNums = (1:length(ims.fibSegs));
        end
        
        curv_all = [ims.fibSegs(:).curv];
        xy_sel = [ims.fibSegs(segNums).xy];
        curv_sel = [ims.fibSegs(segNums).curv];
        [him, hs, hax, hf, hc] = ...
            curv_plot(...
            ims.gray,...
            xy_sel,...
            curv_sel,...
            [min(curv_all),max(curv_all)]);
        
    case 'fib'
        
        if exist('segNums')~=1
            segNums = (1:length(ims.Fibers));
        end
        
        curv_all = [ims.Fibers(:).curv];
        xy_sel = [ims.Fibers(segNums).xy];
        curv_sel = [ims.Fibers(segNums).curv];
        [him, hs, hax, hf, hc] = ...
            curv_plot(...
            ims.gray,...
            xy_sel,...
            curv_sel,...
            [min(curv_all),max(curv_all)]);
end

end