% First have to do Run_Train, sens_analysis, then run GTFiber at the
% minimum parameters and load('sf2debug_2')...

%% FLD histograms for GTFiber and Manual tracing
fld_fa=FiberLengthsFA(fa_path);
figure;histogram(ims.FLD,(0:100:2500))
xlabel('Fiber Length (nm)')
ylabel('# of Fibers')
title('GTFiber')
ax=gca;
ax.FontSize=16;
ax.YLim=[0 50];f=gcf;
f.Position=[716.2000  321.0000  492.8000  418.4000];
figure;histogram(fld_fa,(0:100:2500))
xlabel('Fiber Length (nm)')
ylabel('# of Fibers')
title('Manual Tracing')
ax=gca;
ax.FontSize=16;
ax.YLim=[0 50];
f=gcf;
f.Position=[716.2000  321.0000  492.8000  418.4000];

%% Fiber Vec Plots for each
FiberVecPlot(ims,'gt');
FiberVecPlot(fa_path,'fa');

%% Sfull, DecayLen histograms

[op_fa,fibLenDens]=op2d_FA(imageData)
figure; histogram([res_end(:).Sfull],10)
ax=gca
line([op_fa.Sfull,op_fa.Sfull],ax.YLim,'LineWidth',2,'Color',[0.2 0.2 0.2])
ax.FontSize=16;
xlabel('\it{S}_{full}')
ylabel('# of runs')

figure; histogram([res_end(:).decayLen],10)
ax=gca
line([op_fa.decayLen,op_fa.decayLen],ax.YLim,'LineWidth',2,'Color',[0.2 0.2 0.2])
ax.FontSize=16;
xlabel('Decay Length (nm)')
ylabel('# of runs')

figure; histogram([res_end(:).fibLengthDensity],10)
ax=gca
line([fibLenDens,fibLenDens],ax.YLim,'LineWidth',2,'Color',[0.2 0.2 0.2])
ax.FontSize=16;
xlabel('Fiber Length Density (1/micron)')
ylabel('# of runs')
