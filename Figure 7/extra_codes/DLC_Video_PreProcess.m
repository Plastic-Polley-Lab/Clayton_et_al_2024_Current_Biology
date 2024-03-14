% path = uigetdir('Z:\KeChen\Analysis\','Select the analyzed video folders');
function DLC_Video_PreProcess(path)
cd(path)
file = dir('*.csv');
for i = 1:length(file)
    T = readtable(file(i).name);
    data(i).name=file(i).name;
    data(i).Variables = T([1,2],:);
    temp = load_DLC_csv(file(i).name);
    
    % the x position: 2:3:end; y position: 3:3:end; likelihood position:
    % 4:3:end
    
    for j = 1:size(temp,1) % loops through each frame
        loc(j).x = temp(j, 2:3:end);
        loc(j).y = temp(j, 3:3:end);
        loc(j).likelihood = temp(j, 4:3:end);
        
    end
    
    % separate the pupil and snout
    for j = 1:size(temp, 1)
        data(i).snout(j).x = loc(j).x(end);
        data(i).snout(j).y = loc(j).y(end);
        data(i).snout(j).likelihood = loc(j).likelihood(end);
        data(i).pupil(j).x = loc(j).x(1: end-1);
        data(i).pupil(j).y = loc(j).y(1: end-1);
        data(i).pupil(j).likelihood = loc(j).likelihood(1: end-1);
    end
    clear loc
    
end
%%
for i = 1:length(data)
    fprintf('Processing Video # %d\n', i)
    [distance_ap, pupil] = DLC_pupil(data(i).pupil);
    data(i).pupil_fit = pupil;
    data(i).distance_ap = distance_ap;
    clear distance_ap pupil
end
save('Summary_data', 'data')