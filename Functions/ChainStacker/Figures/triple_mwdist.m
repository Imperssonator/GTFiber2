% function f = triple_mwdist()

navy = [31; 119; 180]./255;
green = [44; 160; 44]./255;
red = [214; 39; 40]./255;

[mu10,sig10] = logn_dist_params(10000,10000*2.25);
[mu25,sig25] = logn_dist_params(25000,25000*2.25);
[mu40,sig40] = logn_dist_params(40000,40000*2.25);

pdf10 = @(x) logn_pdf(x,mu10,sig10);
pdf25 = @(x) logn_pdf(x,mu25,sig25);
pdf40 = @(x) logn_pdf(x,mu40,sig40);

Limits = [1 100000];

f = figure;
ax = gca;
hold(ax,'on')
fplot(pdf10,Limits);
fplot(pdf25,Limits);
fplot(pdf40,Limits);

Lines = ax.Children;
Lines(3).Color = red;
Lines(2).Color = green;
Lines(1).Color = navy;

for i = 1:length(Lines)
    Lines(i).LineWidth = 2;
end

ax.FontSize = 16;
ax.Box='on';
ax.LineWidth = 0.75;
ax.XTickLabels= ...
    cellfun(...
    @(x) num2str(10*str2num(x),3),...
    ax.XTickLabels,...
    'UniformOutput',false);
ax.YTickLabels= ...
    ax.YTickLabels;

% Make Second x-axis for inter-fiber distance
ax2 = copyobj(ax,f);
ax2.XAxisLocation='top';
ax2.XTickLabels= ...
    cellfun(...
    @(x) num2str(1000*0.38*str2num(x)/166,3),...
    ax2.XTickLabels,...
    'UniformOutput',false);