function out = findLastDot(FilePath)

DotInd = regexp(FilePath,'[\.]');
out = DotInd(end);

end