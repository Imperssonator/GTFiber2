function ims = Persist_Len(ims)

% Get uniformly distributed data of fibers in nm
step = ims.settings.fiberStep_nm;
data = {ims.Fibers(:).xy_nm};

procStepNum = max(cellfun(@length, data)) - 1;

xdata = step*(1:procStepNum);
ydata = zeros(1, procStepNum); % mean square end-to-end distance
count = zeros(1, procStepNum); % weight of each graph point

for k = 1:length(data)    
    xy = data{k}; % coordinates of the current fiber
    
    n = length(xy); % number of points in the current fiber
    
    shiftNum = min(procStepNum, n-1);
    
    count(1:shiftNum) = count(1:shiftNum) + (n-1:-1:n-shiftNum);
    
    % Cycle for different separations between points along a fiber
    for l = 1:shiftNum
        ydata(l) = ydata(l) + ...
            sum( sum( ( xy(:,1:end-l) - xy(:,1+l:end) ).^2 ) );
    end
end

ydata = ydata./count;

% Just do the fit for the full dataset. No reason to restrict it.
fit_fun = @(lp,x) 4*lp*(x-2*lp*(1-exp(-x/(2*lp))));
LP = lsqcurvefit(fit_fun,1000,xdata,ydata);

% figure; hold on
% plot(xdata,ydata)
% fplot(@(x) fit_fun(lp(end),x),[0 max(xdata)])

ims.lp=LP;

end


