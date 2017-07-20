dirPath='/Users/Imperssonator/Documents/MATLAB/GTFiber-Mac/Example Images/Wires';
ad=pwd;
cd(dirPath)
D = dir('*.jpg');
cd(ad)

rx = '(?<=_)\d+(?=_res)';   % finds any number of digits between _ and _res, which is the width in inches

Params = [...
    10;
    30;
    2;
    30;
    3000;
    100;
    400;
    0.1;
    3];

for i = 1:length(D)
    disp(D(i).name)
    D(i).width_in = regexp(D(i).name,rx);
    D(i).width_nm = D(i).width_in * 1000;
    D(i).path = [dirPath, '/', D(i).name];
    D(i).results = get_results(Params,D(i).path,D(i).width_nm);
    disp(D(i).results)
end

%% Analysis

for i = 1:length(D);
    D(i).Sfull =...
        D(i).results.op2d.Sfull;
    D(i).decayLen =...
        D(i).results.op2d.decayLen;
    D(i).meanLen =...
        mean(D(i).results.FLD);
    D(i).medLen =...
        median(D(i).results.FLD);
    D(i).fibLengthDensity =...
        D(i).results.fibLengthDensity;
    D(i).fiberCount = ...
        length(D(i).results.FLD);
end

% Fiber count histogram
fiberCounts = [D(:).fiberCount];
figure;histogram(fiberCounts,(min(fiberCounts)-0.5:max(fiberCounts+0.5)))
ax=gca;
ax.XTick=ax.XTick(2:2:end);
ax.FontSize=16;
xlabel('Fiber Count')
ylabel('Number of Test Images')
hold on
hl=plot([length(TrueDist),length(TrueDist)],[0,max(ax.YTick)],'-r');
hl.LineWidth = 2;

fiberLengths = [D(:).meanLen];
figure;histogram(fiberLengths,20);
ax2=gca;
hold on
ax2.FontSize=16;
xlabel('Mean Fiber Length (inches)')
ylabel('Number of Test Images')
hl2=plot([2586,2586],[0,max(ax2.YTick)],'-r');
hl2.LineWidth = 2;
ax2.XTickLabels=cellfun(@(x) num2str(str2num(x)/1000),ax2.XTickLabels,'UniformOutput',false);

load('IMG_1003_12.mat')
figure;
histogram(imageData.length_nm,(0:500:7000));
ax3=gca;
ax3.FontSize=16;
xlabel('Fiber Length (inches)')
ylabel('Number of Fibers')
ax3.XTickLabels=cellfun(@(x) num2str(str2num(x)/1000),ax3.XTickLabels,'UniformOutput',false);
