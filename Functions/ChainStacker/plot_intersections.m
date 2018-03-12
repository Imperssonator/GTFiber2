

% [3,55,64]
% for SonAge example

% [384 385 386] for 'low density' blade example
% [249 251 264] for 'high density' MFU

figure; hold on;
f=249;
plot(ims.Fibers(f).xy_nm(1,:)',ims.Fibers(f).xy_nm(2,:)','-k')
plot(ims.Fibers(f).xing_pts(:,1),ims.Fibers(f).xing_pts(:,2),'or')
ax=gca;
ax.get
axis equal
ax.YDir='reverse';
for f = [251,264];
chains=ims.Chains(ims.ChainLabels==f,:);
for c = 1:length(chains);
plot(ax,chains(c,[1,3]),chains(c,[2,4]),'-b');
end
end

ax.Visible='off'
ax.Position=[0 0 1 1]
ax.Children(end-1).MarkerSize=10
ax.Children(end).LineWidth = 2