function ensure_dir(target_file)

parts=regexp(target_file,filesep,'split');
target_dir = fullfile(parts{1:end-1});

if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

end