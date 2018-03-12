function out = Mn(dist)
% calculate Mn or DPn for a vector of molecular weights or DP's

num = 0;
den = length(dist);

for i = 1:length(dist)
    num = num+dist(i);
end

out = num/den;

end