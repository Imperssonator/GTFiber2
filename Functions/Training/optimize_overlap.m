function [] = optimize_overlap(im_path,fa_path,nmWid,saveFile)

x0 = [...
    10;
    30;
    3;
    30;
    1500;
    80;
    80;
    0.1;
    3];

lb = [...
    1;
    1;
    0.5;
    1;
    100;
    0;
    0;
    0;
    2];

ub = [...
    100;
    100;
    10;
    100;
    5000;
    300;
    500;
    0.2;
    10];

del = [...
    5;
    5;
    1;
    5;
    300;
    10;
    10;
    0.01;
    0.5];

ofun = @(p) objective_overlap(p,im_path,nmWid,fa_path);

opt_params = sensitive_optimization(ofun,x0,lb,ub,del);
ofun(opt_params);

load('Optim_Results')
load('sf2debug_2')

% Results we are interested in studying
fields = {...
    'Sfull',...
    'decayLen',...
    'meanLen',...
    'fibLengthDensity'};

% Results from last iteration of optimization
res_end = sens_res(:,end);
[~,best] = min([res_end(:).out]);

%% Results from FiberApp
fld_fa=FiberLengthsFA(fa_path);
load(fa_path)
[op_fa,fibLenDens]=op2d_FA(imageData);
res_fa = struct( ...
    'Sfull',op_fa.Sfull,...
    'decayLen',op_fa.decayLen,...
    'meanLen',mean(fld_fa),...
    'medLen',median(fld_fa),...
    'fibLengthDensity',fibLenDens);

for i = 1:length(res_end);
    res_end(i).Sfull =...
        res_end(i).results.op2d.Sfull;
    res_end(i).decayLen =...
        res_end(i).results.op2d.decayLen;
    res_end(i).meanLen =...
        mean(res_end(i).results.FLD);
    res_end(i).medLen =...
        median(res_end(i).results.FLD);
    res_end(i).fibLengthDensity =...
        res_end(i).results.fibLengthDensity;
end

% Structure that will contain the partial derivative
% of result j w.r.t. parameter i
sens=struct();
np = size(param_mat,1);   % Number of parameters

for i = 1:np
    for j = 1:length(fields)
        % difference between results normalized by value of result at
        % center divided by difference in parameter normalized by its
        % own value
        
        sens(i).(fields{j}) = ...
            ( ...
            ( res_end(i+np).(fields{j}) ...
            - res_end(i).(fields{j}) ) ...
            / res_end(end).(fields{j}) ...
            ) ...
            / ...
            ( ...
            ( res_end(i+np).results.Params(i) ...
            - res_end(i).results.Params(i) ) ...
            / res_end(end).results.Params(i) ...
            );
    end
end

save(saveFile)

end