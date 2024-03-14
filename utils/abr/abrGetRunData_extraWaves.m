function [runData] = abrGetRunData_extraWaves(runPath)

%ABRGETRUNDATA is used to get frequency, threshold, and wave 1 amplitude
%information from a single file

%%%Input Variables:
%runPath: path to the text file to be parsed

%%%Output Variables:
%runData: structure with all information for this run

%%
%Open the text file
fp = fopen(runPath, 'rt');

%Find the threshold and frequency tested from the first lines
thresh = fscanf(fp, 'Threshold (dB SPL): %f\n');
freq = fscanf(fp, 'Frequency (kHz): %f\n');

%Read through the lines until reaching the table with the peak info
disp(runPath)
lineInfo = 'xyz123';
while ~strcmp(lineInfo(1:5),'Level')
    lineInfo = fgetl(fp);
    lineInfo
end

%Find the correct info from the current line
headings = strsplit(lineInfo,'\t');
idx_level = find(strcmp({headings{:}}, 'Level')); %Level column
idx_amplitude = find(strcmp({headings{:}}, 'P1 Amplitude')); %Wave 1 amplitude column
idx_baseline = find(strcmp({headings{:}}, '0.3msec Avg')); %Baseline for this trace
idx_n1 = find(strcmp({headings{:}}, 'N1 Amplitude')); %For peak to peak amplitude

%Other waves, added 11/4/22
idx_p2 = find(strcmp({headings{:}}, 'P2 Amplitude')); 
idx_n2 = find(strcmp({headings{:}}, 'N2 Amplitude')); 
idx_p3 = find(strcmp({headings{:}}, 'P3 Amplitude')); 
idx_n3 = find(strcmp({headings{:}}, 'N3 Amplitude')); 
idx_p4 = find(strcmp({headings{:}}, 'P4 Amplitude')); 
idx_n4 = find(strcmp({headings{:}}, 'N4 Amplitude')); 
idx_p5 = find(strcmp({headings{:}}, 'P5 Amplitude')); 
idx_n5 = find(strcmp({headings{:}}, 'N5 Amplitude')); 

%Now read through the table section of the document and get appropriate
%info
lineInfo = fgetl(fp);
count = 1;
while lineInfo ~= -1 %Read through the lines until the end
    data = strsplit(lineInfo,'\t');
    levels(count,1) = str2double(data{idx_level}); %Grab the level and amplitude info
    raw_amp(count,1) = str2double(data{idx_amplitude});
    
    %Get different versions of the wave 1 amplitude
    pk2pk_amp(count,1) = raw_amp(count,1) - str2double(data{idx_n1});
    base_amp(count,1) = raw_amp(count,1) - str2double(data{idx_baseline});
    
    %Get other wave amplitudes, added 11/4/22
    w2_amp(count,1) = str2double(data{idx_p2}) - str2double(data{idx_n2});
    w3_amp(count,1) = str2double(data{idx_p3}) - str2double(data{idx_n3});
    w4_amp(count,1) = str2double(data{idx_p4}) - str2double(data{idx_n4});
    w5_amp(count,1) = str2double(data{idx_p5}) - str2double(data{idx_n5});
    
    lineInfo = fgetl(fp); count = count + 1;
end

%Orient vectors correctly
levels = flipud(levels); raw_amp = flipud(raw_amp); pk2pk_amp = flipud(pk2pk_amp); base_amp = flipud(base_amp);
w2_amp = flipud(w2_amp); w3_amp = flipud(w3_amp); %Added 11/4/22
w4_amp = flipud(w4_amp); w5_amp = flipud(w5_amp);

fclose(fp);

%Store the data in a structure
runData = struct;
runData.freq = freq; runData.threshold = thresh; runData.levels = levels;
runData.raw_amp = raw_amp; runData.pk2pk_amp = pk2pk_amp; runData.base_amp = base_amp;
runData.w2_amp = w2_amp; runData.w3_amp = w3_amp; %Added 11/4/22
runData.w4_amp = w4_amp; runData.w5_amp = w5_amp; 
end
