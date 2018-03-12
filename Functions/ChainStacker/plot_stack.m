function [] = plot_stack(stack)

f=figure;
ax=gca;
hold on

for i = 1:length(stack.chains)
    plot(ax,...
        [stack.jog(i)-stack.chains(i)/2,stack.jog(i)+stack.chains(i)/2],...
        [i*0.38,i*0.38],...
        '-b')
end

axis equal

ax.Position = [0 0 1 1];
ax.YLim = [0 500];
ax.Visible='off';

end