function out = Mw(dist)
% calculate Mw or DPw for a vector of molecular weights or DP's

num = 0;
den = 0;

for i = 1:length(dist)
    num = num+dist(i)*dist(i);
    den = den+dist(i);
end

out = num/den;

end