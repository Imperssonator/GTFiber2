function segx = get_seg_intersect(ims,s)

% Next to-do: get threshold values from CEDbw with seek_inds to quantify "gap"

[m,n] = size(ims.gray);

seek = round(distributePoints(ims.fibSegs(s).xy_seek,1));
x=seek(1,:); y=seek(2,:);
badX= (x>n|x<1);
badY= (y>m|y<1);
seek = seek(:,~(badX|badY));  % Nothing outside of the image
seek = flipud(seek)';   % Now i,j
seek_inds = sub2ind(size(ims.gray),seek(:,1),seek(:,2));

% Create black and white image
seek_img = zeros(size(ims.gray));
seek_img(seek_inds)=1;
seek_img = logical(imdilate(seek_img,strel('disk',2)));

seg_intersect=ims.SegLabels(seek_img);
unique_intersect=unique(seg_intersect);
unique_intersect((unique_intersect==s)|(unique_intersect==0))=[];
intersect_fracs = arrayfun(@(x) sum(seg_intersect==x)/length(ims.fibSegs(x).sortPixInds),...
    unique_intersect);
segx = unique_intersect(intersect_fracs>0.7);



end

function sortPixInds = sort_segs(fibSegs,mainSeg)

% First, get all of the endpoints of the matching segments
xy_ends = [];
seg_sort = [];
for i = 1:length(fibSegs);
    xy_ends = [xy_ends;
        fibSegs(i).sortPixSubs(1,:);
        fibSegs(i).sortPixSubs(end,:)];
    seg_sort = [seg_sort;
        [i 1];
        [i 2]];
end
xy_ends=flipud(xy_ends');

% Now, find their coordinates in terms of the arc length of the seeking
% contour
xy_arc = zeros(size(xy_ends,2),1);
xy_seek = fibSegs(1).xy_seek;
for i = 1:length(xy_arc);
    xy_dist = sum((xy_seek-repmat(xy_ends(:,i),1,size(xy_seek,2))).^2,1);
    [~, xy_arc(i)] = min(xy_dist);
end



[ends_sort,orig_ends] = sort(xy_arc);

seg_sort = seg_sort(orig_ends,:);

sortPixInds = [];
for i = 1:2:size(seg_sort,1);
    if seg_sort(i,2)==1;
        sortPixInds = [sortPixInds;
            fibSegs(seg_sort(i,1)).sortPixInds];
    else
        sortPixInds = [sortPixInds;
            flipud(fibSegs(seg_sort(i,1)).sortPixInds)];
    end
end

end
