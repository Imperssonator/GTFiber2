function hf = FLD_hist(ims,figSave)

if exist('figSave','var')~=1
    figSave=0;
end

hf=figure;
ax=gca;
histogram(ims.FLD,30);
ax.FontSize=20;
xlabel('Fiber Length (nm)');
ylabel('Number of Fibers');

flfont=16;
flpos=[0.55, 0.87];
edgedark = 0;
edgewidth = 0.75;

htex = text('Units', 'normalized', 'Position', flpos, ...
    'BackgroundColor', [1 1 1], ...
    'String', {['<L>= ', num2str(mean(ims.FLD)), ' nm'],...
               ['std. dev.= ', num2str(std(ims.FLD)), ' nm'],...
               ['{\it\rho}_{FL}= ' num2str(ims.fibLengthDensity), ' µm^{-1}']}, ...
    'FontSize', flfont,...
    'EdgeColor', edgedark*[1 1 1],...
    'LineWidth', edgewidth);


if figSave
    F = getframe(hf);
    Fim = F.cdata;
%     Fres = imresize(Fim,[h, w]);
    imwrite(Fim, [ims.figSavePath, '_FLD.tif']);
    close(hf)
end

end