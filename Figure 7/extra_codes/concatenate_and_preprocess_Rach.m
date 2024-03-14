function concatenate_and_preprocess_Rach(analdir,files,penetration_number, rawDataPrefix,NChan,chans,region,blockNums, raw_waveform_path)
%function that takes in a cell array of files and sorts them together

disp('Preparing to sort bulk data...');
disp('Files to be sorted are:');
for i=1:length(files)
    disp(sprintf('...%s',files{i}));
end

%this analysis directory will never change
% raw_analdir = '/Volumes/RSW_Data/Raw_Data_mat/';
% raw_analdir = [analdir '/Raw_Waveform/'];
raw_analdir = [raw_waveform_path,'\'];

filename = [rawDataPrefix '_complete'];

nchan=NChan;       %number of channels 
sr=24414.0625;    %sampling rate (should always be 24414.06 so as not to drop any samples)

blocks=[];
tanks=[];
save_files = [];
count=1;
for i=1:length(files)
    file = files{i};
    load(file)
    
%     s = discrimSettings.BlockName
%     a = sscanf(s,'D:\Data\Tanks\MMA%d\MMA%d-%d-%d')% change based on what's stored in the impale file discrimSettings 
    %catch if is a Darwin file
%     if ~isempty(strfind(file,'Darwin'))
%         for j=1:length(GenInfo)
%             s_ = GenInfo(j).Block;
%             a_ = sscanf(s_,'D:\Data\Tanks\MMA%06d\MMA%06d-%d-%d');
%             blocks(count) = a_(4);
%             penetrations(count) = a_(3);
%             tanks(count) = a_(1);
%             save_files{count} = files{i};
%             count=count+1;
%         end     
%     else    
        blocks(count) = i;
        penetrations(count) = str2num(rawDataPrefix(end));
%         tanks(count) = str2num(rawDataPrefix(end-7:end-2));
        save_files{count} = files{i};
        count=count+1;      
%     end  
end


%NOW DO THE SORTING
filepath = [analdir sprintf('/Sorted_Data/%s/Penetration_%d/',region,penetration_number)];
mkdir(filepath)


%first load in the data
disp('Loading data...');

data = [];
store_ticks = [];
timepoints = [];
tankname = split(analdir, '\');
tankname = char(tankname(end));
for i=1:length(blockNums)
    waveform_channels = [];
    disp(sprintf('...block %d/%d',i,length(blockNums)));
    for j=1:nchan
        disp(sprintf('...... channel %d',chans(j)));
        load([raw_analdir sprintf('waveform_tank_%s_pen_%d_block_%d_channel_%d',tankname,penetrations(i),blockNums(i),chans(j))]);
        waveform_channels(j,:) = results.waveform .* 1e6; 
        %multiply by 1e6 if data saved as float (to get into uv)
        %divide by 4 if data saved as int16 (to get into uv)
    end
    
    waveform_channels = int16(waveform_channels);
    
    %concat one channel, just to get lengths correct
    curr_time = length(data)./sr;
    data = [data waveform_channels(1,:)];
    
    disp(sprintf('.........starting to write block %d to disk',i))
    if i==1
        fp = fopen([filepath 'data.dat'], 'wb');
        fwrite(fp, waveform_channels(:), 'int16');
        fclose(fp);
    else
        fp = fopen([filepath 'data.dat'], 'ab');
        fwrite(fp, waveform_channels(:), 'int16');
        fclose(fp);    
    end
    disp(sprintf('.........finished writing block %d to disk',i))

    load([raw_analdir sprintf('data_tank_%s_pen_%d_block_%d.mat',tankname,penetrations(i),blockNums(i))]);
    grab_ticks = results.ticks + curr_time;
    store_ticks = [store_ticks grab_ticks];
    timepoints = [timepoints curr_time];
end

% add common average to filter out artifact by Ke
cd(filepath)
applyCARtoDat('data.dat', nchan);

savename = [filename '_timepoints'];
save([filepath savename],'timepoints');

savename = [filename '_ticks'];
save([filepath savename],'store_ticks');

savename = [filename '_files'];
save([filepath savename],'save_files');

% fp = fopen([filepath 'data.dat'], 'wb');
% fwrite(fp, data(:), 'int16');
% fclose(fp);


clear data waveform_channels %for memory










