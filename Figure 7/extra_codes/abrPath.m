function [finalPath] = abrPath(mouse, session, subFolder)

%ABRPATH is used to get a path to a mouse's abr folder

%%%Input Variables:
%mouse: mouse name
%session: session number
%subFolder: which sub-folder within the 'ABR Data' folder that the data is found in

%%%Output Variables:
%finalPath: path to the specified file of interest


%%
%First define the overall directory

if nargin < 2
    startPath = 'D:\ABR Data';
else
    startPath = strcat('D:\ABR Data\',subFolder);
end


%Now use the mouse, session, and run to find the proper folder
format = '\\%s\\%s%d\\';
endPath = sprintf(format,strcat(mouse,' DATA'),mouse,session);

%Final path
finalPath = strcat(startPath,endPath);


end
