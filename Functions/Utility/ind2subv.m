function v = ind2subv(siz,ndx) 
[out{1:length(siz)}] = ind2sub(siz,ndx); 
v = cell2mat(out);