function out = findLastSlash(FilePath)

if ispc
    SlashInd = regexp(FilePath,'[\\]');
else
    SlashInd = regexp(FilePath,'[\\/]');
end

if isempty(SlashInd)
    out = 0;
else
    out = SlashInd(end);
end

end