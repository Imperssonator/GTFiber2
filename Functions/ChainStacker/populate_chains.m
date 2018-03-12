function [ims, Fim] = populate_chains(ims,Mn,PDI,figSwitch)

% Populate Chains
%
% Inputs: ims structure with Fibers field
%
% Outputs: ims structure with new field "Chains", which is a 2x2xN matrix
% of all of the endpoints of individual polymer chains in the fibers in the
% film.

% The format of each entry in the Chains matrix is like this:
% | x1  x2 |
% | y1  y2 |


% Set up figure if plotting
if exist('figSwitch','var')~=1
    figSwitch=0;
end

if figSwitch
    f1 = figure('Visible','on');
    f1.Position = [100 100 800 800];
    f1.Color = [1 1 1];
    % f1.InvertHardcopy='off';
    ax = axes('parent',f1);
    hold(ax,'on');
    ax.YDir='reverse';
    ax.XLim = [0, size(ims.img,2)*ims.nmPix];
    ax.YLim = [0, size(ims.img,1)*ims.nmPix];
    ax.Visible = 'off';
    % ax.PlotBoxAspectRatio = [1 1 1];
    ax.Position = [0 0 1 1];
end


% Set up molecular weight distribution for sampling
pi_d = 0.38; mon_weight=166; extra_length=0;
[mu,sig] = logn_dist_params(Mn,Mn*PDI);


% Initialize storage variables
Chains = [];
ChainLabels = [];


for f = 1:length(ims.Fibers)
    disp(f)
    % Each fiber is composed of discrete straight segments defined by their
    % endpoints in x,y space. For each segment, we will sample the
    % molecular weight distrubtion defined above by Mw and Mn, and add
    % chains of the sampled lengths at intervals along the fiber segments
    % spaced by the pi-pi stacking distance such that they lie
    % perpendicular to the orientation of that fiber segment, and are
    % centered on the backbone of that fiber segment.
    
    for seg = 1:(size(ims.Fibers(f).xy,2)-1)
        
        xy_nm = ims.Fibers(f).xy_nm; % .* ims.nmPix;
        
        pt1 = xy_nm(:,seg);
        pt2 = xy_nm(:,seg+1);
        
        stack_vec = pt2-pt1;
        stack_pts = (extra_length:pi_d:norm(stack_vec))';
        
        % We don't want the pi-pi stacking distance to reset for every segment -
        % extra_length is how much to offset the start of the next segment's stack
        % (in nm)
        extra_length = pi_d-(norm(stack_vec)-stack_pts(end));
        
        stack_fracs = stack_pts/norm(stack_vec);
        stack_pts_2d = repmat( pt1', size(stack_pts,1), 1 ) + ...
                       repmat( stack_vec', size(stack_pts,1), 1 ) ...
                       .* repmat( stack_fracs,1,2 );
        
        chain_lengths = ceil( ...
                             lognrnd( mu, sig, size(stack_pts,1), 1 ) ...
                             ./ mon_weight ...
                            ) ...
                            .*pi_d;
        half_chain_lengths = chain_lengths./2;
        
        % Here we are adding each half-chain-length perpendicularly outward
        % from each stack_pt
        
        ortho_vec = [-stack_vec(2); stack_vec(1)];
        ortho_unit_vec = ortho_vec./norm(ortho_vec);
        
        % The stack will be an Nx4 stack of chain endpoints:
        % | x1_end1 y1_end1 x1_end2 y1_end2 |
        % | xn_end1 yn_end1 xn_end2 yn_end2 |
        
        stack = zeros(size(stack_pts,1),4);
        stack(:,1) = stack_pts_2d(:,1)-half_chain_lengths.*ortho_unit_vec(1);
        stack(:,2) = stack_pts_2d(:,2)-half_chain_lengths.*ortho_unit_vec(2);
        stack(:,3) = stack_pts_2d(:,1)+half_chain_lengths.*ortho_unit_vec(1);
        stack(:,4) = stack_pts_2d(:,2)+half_chain_lengths.*ortho_unit_vec(2);
        
        StackLabels = ones(size(stack,1),1).*f;
        ChainLabels = cat(1,ChainLabels,StackLabels);
        Chains = cat(1,Chains,stack);
        
        if figSwitch
            if 1 %ismember(f,[1, 2, 3, 4, 5, 6, 7, 8, 9])
                for i = 1:size(stack,1)
                    plot(ax,stack(i,[1,3]),stack(i,[2,4]),'-b');
                end
            end
        end
        
    end
end

ims.Chains = Chains;
ims.ChainLabels = ChainLabels;

if figSwitch
    F = getframe(f1);
    Fim = F.cdata;
end

end