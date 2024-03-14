function [allData] = abrGetSessionData(finalPath)

%ABRGETSESSIONBDATA is used to get all information from one session

%%%Input Variables:
%path: path to the session folder

%%%Output Variables:
%allData: structure array of information for all frequencies tested

%%

%Go to the correct directory and grab those files
files_in_dir = dir(finalPath);

%Find all files that indicate they are analyzed
toAnalyze = '-analyzed'; count = 1;

for i = 1:length(files_in_dir) %Loop through all files in directory
    
    currFile = files_in_dir(i).name %Grab file name
    
    if contains(currFile, toAnalyze) %Only look at analyzed ABR files
        
        %Grab the path to this file
        runPath = strcat(finalPath,currFile);
        %Get the run data
        allData(count) = abrGetRunData_extraWaves(runPath);
        
        count = count + 1;
    end
end
        