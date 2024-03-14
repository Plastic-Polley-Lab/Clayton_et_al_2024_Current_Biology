% get all sev files
clear
files = dir('*.sev')

% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
    parts = split(f,'_')
    f1 = f(1:12)
    f2 = char(parts(5))
    f3 = char(parts(6))
    run_num = 10;
    % Convert to number
    % If numeric, rename
    newnames = [f1, '_', f1, '-2-', num2str(run_num), '_', f2, '_', f3, '.sev']
    movefile(files(id).name, newnames);
end
%%