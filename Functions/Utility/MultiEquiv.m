function out = MultiEquiv(M,V)

%Multiple Equivalence Test
% Given a vector of values, V, return out, a matrix the same size as M,
% with ones where a value of M is equal to ANY of the values in V.

out = zeros(size(M));

for i = 1:length(V)
    
    out = out+double(M==V(i));
    
end

out = logical(out);

end