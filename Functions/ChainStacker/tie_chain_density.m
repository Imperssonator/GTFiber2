function ims = tie_chain_density(ims)

% Get tie chain length density for a structure in which chains have been
% populated

% Hard coded:
% edge_buffer: fraction of image width used as buffer zone to exclude
% calculation of tie chain density for fibers subject to edge effects
edge_buffer=0.05;

for f = 1:length(ims.Fibers)
    
    disp(f)
    
    % Do some pre-processing stuff that probably should've been done
    % earlier
    xy_nm = ims.Fibers(f).xy .* ims.nmPix;
    ims.Fibers(f).xy_nm = xy_nm;
    ims.Fibers(f).length_nm = sum(sqrt(sum(diff(ims.Fibers(f).xy_nm,1,2).^2,1)));

    % Check if any part of fiber lies within edge buffer
    if any(xy_nm(:)<(edge_buffer*ims.nmWid)) ||...
            any(xy_nm(:)>((1-edge_buffer)*ims.nmWid))
        disp('too close to edge')
        ims.Fibers(f).xing_pts=[];
        ims.Fibers(f).tieChainDensity=[];
        continue
    end
    
    % Get 'bounding box' coordinates of each fiber to expedite chain
    % screening
    fmax_x = max(xy_nm(1,:));
    fmin_x = min(xy_nm(1,:));
    fmax_y = max(xy_nm(2,:));
    fmin_y = min(xy_nm(2,:));
    
    % Initialize matrix to store crossing points of tie chains for each
    % fiber
    xing_pts = [];
    
    % Exclude chains in current fiber
    chains_other = ims.Chains(ims.ChainLabels~=f,:);
    chains_out = (chains_other(:,1)<fmin_x & chains_other(:,3)<fmin_x) |...
                 (chains_other(:,1)>fmax_x & chains_other(:,3)>fmax_x) |...
                 (chains_other(:,2)<fmin_y & chains_other(:,4)<fmin_y) |...
                 (chains_other(:,2)>fmax_y & chains_other(:,4)>fmax_y);
    chains_other = chains_other(~chains_out,:);

    xy_nm_t = xy_nm';
    xy_nm_linesegs = [xy_nm_t(1:end-1,:), xy_nm_t(2:end,:)];
    inter_struct = lineSegmentIntersect(xy_nm_linesegs,chains_other);
    xing_pts = [inter_struct.intMatrixX(inter_struct.intMatrixX~=0),...
                inter_struct.intMatrixY(inter_struct.intMatrixY~=0)];

    ims.Fibers(f).xing_pts = xing_pts;
    ims.Fibers(f).tieChainDensity = size(ims.Fibers(f).xing_pts,1)/ims.Fibers(f).length_nm;

end

% Compile and store tie chain density distribution as field in ims
in_bounds=arrayfun(@(x) ~isempty(x.tieChainDensity),ims.Fibers);
ims.tc_dist=[ims.Fibers(in_bounds).tieChainDensity]';
    
end