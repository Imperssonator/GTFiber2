function stack = pi_stack(num_chains,Mn,Mw,af_start,len_start,prob_step)

% Grow a pi-stack with both step and condensing growth from a specified
% number of chains with a specified fraction involved in nuclei of a
% specified length

chains_start = round(num_chains * (1-af_start));
nuclei = round(num_chains * af_start / len_start);

% Stack is a structure array that keeps track of every unique unit
% Fields:
%   Chains: vector of chain lengths
%   Jog: vector of center points of each chain relative to center of first
%   chain in unit

[mu,sig] = logn_dist_params(Mn,Mw);
chains = ceil(lognrnd(mu,sig,num_chains,1));
chain_lens = chains./166.*0.38; % Chain lengths in nm

stack = struct('chains',[],'jog',[]);

for i = 1:nuclei
    stack(i).chains = chain_lens((i-1)*len_start+1:i*len_start);
    stack(i).jog = zeros(length(stack(i).chains),1);
end

count=nuclei;
for i = nuclei*len_start+1:length(chains)
    count = count+1;
    stack(count).chains = chain_lens(i);
    stack(count).jog = 0;
end

disp('stack em up!')

while length(stack)>1
    disp(length(stack))
    col_i = randi(length(stack));
    col_j = randi(length(stack));
    while col_j==col_i
        col_j = randi(length(stack));
    end
    
    stack = collide_ij(stack,col_i,col_j,prob_step);
end

plot_stack(stack);

end

function stack = collide_ij(stack,i,j,prob_step)

% reject nucleation events
if length(stack(i).chains)<2 && length(stack(j).chains)<2
    return
% reject step growth with probability prob_step
elseif length(stack(i).chains)<2 || length(stack(j).chains)<2
    if rand>prob_step
        return
    end
end

% Pick a collision point according to a gaussian distribution with std. dev
% = half of the last chain width of i / 4
half_width_end = stack(i).chains(end);
collide_pt = randn()*half_width_end/40;
if collide_pt>half_width_end
    collide_pt=half_width_end;
elseif collide_pt<-half_width_end
    collide_pt=-half_width_end;
end

new_unit = stack(j);
new_unit.jog = new_unit.jog+collide_pt+stack(i).jog(end);

stack(i).chains = [stack(i).chains; new_unit.chains];
stack(i).jog = [stack(i).jog; new_unit.jog];

stack(j) = [];

end