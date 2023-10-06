function newFilename = generate_numbered_filename(baseName, directory)
    existingFiles = dir(fullfile(directory, [baseName, '*'])); % List existing files
    
    [~, name, ext] = fileparts(baseName); % Extract filename parts
    numbering = 1; % Initialize numbering
    
    % Check for clashes and update numbering
    while ~isempty(existingFiles)
        baseName = [name, '_', num2str(numbering), ext];
        existingFiles = dir(fullfile(directory, [baseName, '*']));
        numbering = numbering + 1;
    end
    
    newFilename = fullfile(directory, baseName);
end
