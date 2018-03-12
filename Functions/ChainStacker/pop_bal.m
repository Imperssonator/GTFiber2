function [] = pop_bal(af_start,af_end,chains_total,len_start,len_max,kn,kg1,kg2)

% Population balance simulation

% Start with total_chains number of chains, and a number of nuclei equal to
% (total_chains * start_af) / start_len
% From there, calculate the growth rate of each length bin as
% dnl_(i+j)/dt = kg * [nl_i] * [nl_j]
% But instead of calculating the growth rate of the resulting fiber, we
% will calculate the collision rate of all the segments, and roll for a
% collision, instead of rolling for the growth of a particular bin

% In other words, we're not going to be integrating this as a diff eq. It's
% just going to be a kinetic monte carlo simulation of fiber growth.

len_max = round(len_max);
pop = zeros(len_max,1);

% Start with specified nuclei

pop(1) = chains_total * (1-af_start);
pop(len_start) = round(chains_total * af_start / len_start);

% Fix chain balance
pop = fix_chain_balance(pop,chains_total);

% Calculate current aggregate fraction - simulation runs until
% af_cur==af_end
af_cur = 1-pop(1)/chains_total;

iter = 1;


% Build lookup table for linear index of a low triangle to its subscript
% indices
r2ij = [];
ss = size(pop,1);
for i = 1:ss
r2ij = [r2ij; [(i:ss)',ones(ss-(i-1),1).*i]];
end

% Build the k matrix
k_mat = ones(length(pop),length(pop)) .* kg2;   % Most entries are kg2
k_mat(1,:) = kg1; k_mat(:,1) = kg1;
k_mat(1,1) = kn;    % Probably not going to allow nucleation but that's ok

while af_cur<af_end
    
    iter = iter+1;
    
    % Calculate rates
    % First calculate all of the concentration products:
    C_prod = pop * pop';
    
    Rates = k_mat .* C_prod;
    
    





end


function pop = fix_chain_balance(pop,chains_total)

in_fibers = sum( (2:length(pop))' .* pop(2:end) );

pop(1) = chains_total - in_fibers;

end