for s = 1:length(ims.fibSegs);
    disp(s)
    fiberIntensity = mean(double(ims.gray(ims.segsInit)));
    xy_init=ims.fibSegs(s).xy; %flipud(ims.fibSegs(s).sortPixSubs');
    %     k2_list = 100;
    %     fl_list = zeros(size(k2_list));
    %     for i = 1:length(k2_list);
    params = [0, 1000, 40, 40, 140, fiberIntensity, 10000];
    [xy,curv] = FitSeg(ims,xy_init,ims.gray,ims.settings.fiberStep,params);
    %         fl_list(i) = length(xy);
    %     end
    ims.fibSegs(s).xy_seek=xy;
end

L = zeros(length(ims.fibSegs));
for s = 1:length(ims.fibSegs)
    disp(s)
    true_matches = get_seg_intersect(ims,s);
    L(s,true_matches) = 1;
end

%%

if isfield(ims,'FibersNew')
    ims=rmfield(ims,'FibersNew');
end
% 
% MM = L.*L';
% MM1 = MM*MM;
% MM2 = MM1*MM1;
% MM3 = MM2*MM2;

Li = L+eye(size(L));
Mutuals = logical(L.*L'+eye(size(L)));
[C,IA,IC] = unique(Mutuals,'rows');
MU=double(C)*double(C');

GMU = zeros(size(MU,1),1);
for u = 1:length(GMU);
GMU(u) = all(all( ...
Li(MultiEquiv(IC,find(MU(u,:))),MultiEquiv(IC,find(MU(u,:)))) | ...
Li(MultiEquiv(IC,find(MU(u,:))),MultiEquiv(IC,find(MU(u,:))))' ...
));
end
GMU=logical(GMU);
GC=MU(GMU,GMU);
good2mu = find(GMU);

M_copy = Mutuals;
Mutuals = MU;

stitched = [];
fiber_count = 0;
FiberLabels=zeros(size(ims.SegLabels));
seg2fib = zeros(size(Mutuals,1),1);
for i = 1:size(Mutuals,1)
    segNums = find(MultiEquiv(IC,find(Mutuals(i,:))))';
%     segNums = find(Mutuals(i,:));
    if ~all(ismember(segNums,stitched))
        fiber_count = fiber_count+1;
        fibSegs = ims.fibSegs(segNums);
        [xy,curv] = FitFiber2(ims,fibSegs);
        ims.FibersNew(fiber_count).xy = xy;
        ims.FibersNew(fiber_count).curv = curv;
        ims.FibersNew(fiber_count).fibSegs = segNums;
%         ims.FibersNew(fiber_count).initSeg = i;
        ims.FibersNew(fiber_count).initSeg = i;
        stitched = [stitched, segNums];
        FiberLabels = FiberLabels + double(MultiEquiv(ims.SegLabels,segNums)).*fiber_count;
        seg2fib(segNums) = fiber_count;
    end
end
%
% fiber_count = 1;
% while fiber_count<=length(ims.FibersNew)
%     disp(fiber_count)
%     disp(length(ims.FibersNew))
%     all_others = (1:length(ims.FibersNew));
%     all_others(all_others==fiber_count)=[];
%     all_others_img = ContourImageFiber(ims,all_others,5);
%     this_fiber_img = ContourImageFiber(ims,fiber_count,5);
%     overlap_img = this_fiber_img.*all_others_img;
%     if sum(overlap_img(:))/sum(this_fiber_img(:))>0.7
%         ims.FibersNew(fiber_count) = [];
%     else
%         fiber_count = fiber_count + 1;
%     end
% end
%
FiberColors(ims);

%
% VV=GCAFG(L,70);
%
% for c = 1:max(VV)
%     community = find(VV==c);
%     comm_L = L(community,community);
%     Mutuals = comm_L.*comm_L'+eye(length(community));
% end