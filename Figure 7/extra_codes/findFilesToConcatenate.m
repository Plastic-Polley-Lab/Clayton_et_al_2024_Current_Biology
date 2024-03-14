function fileList= findFilesToConcatenate(rawDataPath, rawDataPrefix)
    % is the raw data path pointing to a directory full of mat files?
	assert(isdir(rawDataPath), 'Data path not found!');
    % make sure the prefix has a wildcard on it to make dir2 happy    
    if(rawDataPrefix(end) ~= '*')
        rawDataPrefix = [rawDataPrefix,'*'];
    end
    % find all mat files that match and sort them lexicographically
    fileList = sort_nat(dir2('/s',rawDataPath, rawDataPrefix,'.mat'));
%     fileList = sort_nat(dir2(rawDataPath, rawDataPrefix,'.mat'));
    assert(length(fileList)>1,'Less than two files found to concatenate. That makes no sense.');
 end