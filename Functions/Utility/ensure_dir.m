function ensure_dir(target_file)

parts=regexp(target_file,filesep,'split');
target_dir = fullfile(parts{1:end-1});
if target_dir(1:5)=='Users'
    % The hackiest hack that ever did hack
    target_dir = [filesep, target_dir];
end
% disp(target_dir)

if ~exist(target_dir, 'dir')
    mkdir(target_dir);
end

end