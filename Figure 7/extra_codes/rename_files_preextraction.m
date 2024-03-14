function rename_files_preextraction(folderpath)
%To remove things like (2) if needed.....

% clearvars;
% folderpath = 'F:\EPHYS\Rach_Tanks\MMA060719-Session1\MMA060719-Session1-Run22';

cd(folderpath);
Files=dir('*.*');
for k=1:length(Files)
   FileNames = Files(k).name;
   splitstr = split(FileNames,' ');
   if (length(splitstr)==2)
       splitstr2 = split(splitstr{2},'.');
       temp = splitstr2{1}(4:end);
       FileNameNew = strcat(splitstr{1},temp,'.',splitstr2{2});
      movefile(fullfile(folderpath,FileNames),fullfile(folderpath,FileNameNew));
   end

end

end