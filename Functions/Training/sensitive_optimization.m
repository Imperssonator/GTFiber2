function opt_params = sensitive_optimization(ofun,x0,lb,ub,del)

max_iter = 10;
cur_iter = 0;
num_params = length(x0);

% Initialize parameter matrix
params = zeros(length(x0),3);
params(:,2) = x0;
params(:,1) = x0-del; params(:,3) = x0+del;
too_low = params(:,1)<lb;
too_high = params(:,3)>ub;
params(too_low)=lb(too_low);
params(too_low,1) = lb(too_low);
params(too_high,3) = ub(too_high);
run = 0;

sens_res = struct();

while cur_iter<max_iter
    cur_iter=cur_iter+1;
    disp(cur_iter)
    
    % This is skipped on the first iter, for some reason I have it
    % happening at the start of the loop though. Probably to avoid
    % instantiating stuff outside the loop. Looks like I'm not even using
    % the regression in any meaningful way, just taking the 'best' (lowest
    % obj. fun. value) setting for each parameter and making that the new
    % center. This algorithm doesn't even have any kind of guarantee of
    % convergence but with this many parameters and such a complex
    % response, this is the only way to do it.
    
    if cur_iter>1
        for p = 1:num_params
            [~,min_pt] = min(param_reg(p).reg.y);
            params_new(p,2)=params(p,min_pt);
        end
        params_new(:,1) = params_new(:,2)-del;
        params_new(:,3) = params_new(:,2)+del;
        too_low = params_new(:,1)<lb;
        too_high = params_new(:,3)>ub;
        params_new(too_low)=lb(too_low);
        params_new(too_low,1) = lb(too_low);
        params_new(too_high,3) = ub(too_high);
        params=params_new;
    end
    
    disp(params)
    
    % This builds a matrix numParams x numParams*2+1, where the last column
    % is the original parameter set, the first 1:numParams columns have
    % each parameter set at its lower value, and the next
    % numParams+1:2*numParams columns have each parameter set at its higher
    % value
    param_mat = repmat(params(:,2),1,num_params*2+1);
    for i = 1:num_params
        param_mat(i,i) = params(i,1);
    end
    for i = 1:num_params
        param_mat(i,i+num_params) = params(i,3);
    end
    
    % sens_res is a structure array that stores the results for each of the
    % numParams*2+1 parameter sets in param_mat
    % sens_res.out has the objective function value, while sens_res.Fibers
    % has the final set of Fibers
    for i = 1:size(param_mat,2)
        run = run+1;
        disp('Run number')
        disp(run)
        disp('----')
        disp('Params')
        disp(param_mat(:,i))
        [sens_res(i,cur_iter).out,sens_res(i,cur_iter).results] = ofun(param_mat(:,i));
    end
    
    % param_reg is a structure array that performs a linear regression on
    % the objective function values for the three levels of each varied
    % parameter
    param_reg = struct();
    for i = 1:num_params;
        param_reg(i).reg = ...
            MultiPolyRegress(...
            params(i,:)',...
            [sens_res(i,cur_iter).out;
            sens_res(num_params*2+1,cur_iter).out;
            sens_res(i+num_params,cur_iter).out],...
            1);
    end
    for i = 1:num_params;
        cv(i) = param_reg(i).reg.CVMAE;
    end
end

[~,best_set] = min([sens_res(:,cur_iter).out]);
opt_params=param_mat(:,best_set);

save(['Optim_Results.mat'],'sens_res','param_mat')

end
