function handles = imshowGT(img,handles,ax_tag)

% hax is the name of a field of handles that may or may not exist that
% references a set of axes which may or may not exist. Check if it exists,
% if not, create it, then put the image on it.

if isfield(handles,ax_tag)
    if isvalid(handles.(ax_tag))
%         disp(1)
%         handles.(ax_tag).Children
        delete(handles.(ax_tag).Children)
        imshow(img,'parent',handles.(ax_tag))
        axis(handles.(ax_tag),'equal');
        handles.(ax_tag).Position = [0 0 1 1];
    else
%         disp(2)
        figure;
        handles.(ax_tag) = axes();
        imshow(img,'parent',handles.(ax_tag))
        axis(handles.(ax_tag),'equal');
        handles.(ax_tag).Position = [0 0 1 1];
    end
else
%     disp(3)
    figure;
    handles.(ax_tag) = axes();
    imshow(img,'parent',handles.(ax_tag))
    axis(handles.(ax_tag),'equal');
    handles.(ax_tag).Position = [0 0 1 1];
end

end