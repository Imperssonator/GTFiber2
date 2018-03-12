function stats = test_dist(DPavg,SD)

sample = lognrnd(DPavg,SD,100000,1);
%disp(length(sample))

% table = tabulate(sample);                   % [DP, count, mol frac]
% out = table(table(:,2)~=0 & table(:,1)>0,:);
% 
DPn = Mn(sample);
DPw = Mw(sample);
PDI = DPw/DPn;

stats = [DPn,DPw];

end
