function handles = plotGT(hax,handles,ax_tag)

% hax is an axis where something has been plotted, but we want it to be on
% one of the GTFiber axes. First we have to clear the existing content,
% then copy the content from hax to the parent of ax_tag

if isfield(handles,ax_tag)
    if isvalid(handles.(ax_tag))
%         disp(1)
        fp = handles.(ax_tag).Parent;
        delete(fp.Children);
        handles.(ax_tag) = copyobj(hax,fp);
    else
%         disp(2)
        fp = figure;
        handles.(ax_tag) = copyobj(hax,fp);
    end
else
%     disp(3)
    fp = figure;
    handles.(ax_tag) = copyobj(hax,fp);
end

end