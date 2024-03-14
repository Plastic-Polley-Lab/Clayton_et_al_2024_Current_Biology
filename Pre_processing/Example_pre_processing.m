%% Baseline Individual Plots
%this script is to be used with raw data. enter which sessions and runs
%correspond to the different time points. If absolute motion energy has not
%been extracted, user will be prompted to draw and roi and wait for
%orofacial energy be calculated for all videos. To use, run full script
%after appropiate changes.

%the abs energy matrix, intensity level and frequency for each trial, and
%startle amplitude will be saved within the processed data folder for each
%mouse to be used in summary plots. Will plot abs motion energy for both
%freqs. Will also plot startle (change in body of script)

%summary plot script: baseline_sum.m (\\apollo\Polley_Lab\Mouse Videography
%Audiogram\Code)
%%

%User change this part:
types = {'pre1';'pre2';'post2h';'post5d';'post2w'};
for blt = 1:length(types)
    baseline_type = types{blt}; %pre2 2hpost 5dpost 2wpost
    mice = [82 84 85];
    
    if strcmp(baseline_type, 'pre1') == 1
        sessions = [5 5 5];
        run = [1 1 1];
    elseif strcmp(baseline_type, 'pre2') == 1
        sessions = [6 6 6];
        run = [1 1 1];
    elseif strcmp(baseline_type, 'post2h') == 1
        sessions = [7 7 7];
        run = [1 1 1];
    elseif strcmp(baseline_type, 'post5d') == 1
        sessions = [8 8 8];
        run = [1 1 1];
    elseif strcmp(baseline_type, 'post2w') == 1
        sessions = [9 9 9];
        run = [1 1 2];
    end
    
    for n = 1:length(mice)
        clearvars -except n mice sessions run baseline_type blt types
        
        addpath(genpath('C:\MouseVideoAudiogram\Code\KC_pilot_code\utils'))
        addpath(genpath('C:\MouseVideoAudiogram\Code\KC_pilot_code'))
        pseudo_cond = 0;%1 or 0 depending on whether it was a pesudo-conditioning session
        
        mouse_num = mice(n);
        session_num = sessions(n);
        run_num = run(n);
        
        %Specify directories to use
        data_root = '//apollo/Polley_Lab/Mouse Videography Audiogram/Data';
        animalID = sprintf('Chat%02d',mouse_num);
        data_dir = fullfile(data_root,sprintf('Chat%02d/Session %d',mouse_num,session_num)); %location of the data (videos)
        aviFolder = fullfile(data_root,sprintf('Chat%02d/Session %d',mouse_num,session_num)); %location of videos
        savingDir = fullfile(data_root,sprintf('Chat%02d/proc_orofacial/Session %d/Run%d',mouse_num,session_num,run_num)); %saving location for processed videos
        
        mkdir(savingDir);
        cd(savingDir);
        
        try
            load(sprintf('Chat%02d_Orofacial.mat',mouse_num));
        catch
            get_orofacial_ROIs_compute_Motion_ChatMice_startle(mouse_num, session_num, run_num)
            cd(savingDir);
            load(sprintf('Chat%02d_Orofacial.mat',mouse_num));
        end
        
        %% Read in tosca trial data
        tic
        try
            %Try loading the variable, will obviously work only if it already exists
            load(fullfile(savingDir,'tl.mat'));
        catch error_info
            fprintf('tl did not exist, creating and saving...');
            %     [d,p] = tosca_read_run(fullfile(data_dir,sprintf('%sAA%d-Session%d-Run%d.txt',behav_vs_phys,mouse_num,session_num, run_num)));
            %     [d,p] = tosca_read_run(fullfile(data_dir,sprintf('AA%s%d-Session%d-Run%d.txt',behav_vs_phys,mouse_num,session_num, run_num)));
            [d,p] = tosca_read_run(fullfile(data_dir,sprintf('Chat%2d-Session%d-Run%d.txt',mouse_num,session_num, run_num)));
            nan_i = [];
            for i = 1:length(d)
                if isnan(d{i}.N)
                    nan_i = [nan_i i];
                end
            end
            for i = 1:length(nan_i)
                d(nan_i) = [];
            end
            %     d = d(1:10);
            tl = tosca_create_log(fullfile(data_dir,sprintf('Chat%2d-Session%d-Run%d.txt',mouse_num,session_num, run_num)),...
                'aviFolder', aviFolder);
            
            % Add analog data read-in
            tl = read_analog_data (tl,d,p);
            
            
            %remove error trials if any, or aborted trials if any
            error_abort_trials = [];
            for trial_num = 1:length(tl.trials)
                if (strcmp(tl.trials{trial_num}.Result, 'Error')||strcmp(tl.trials{trial_num}.Result(end-4 : end), 'Error'))
                    error_abort_trials = [error_abort_trials trial_num];
                end
            end
            for trial_num = 1:length(tl.trials)
                history = tl.trials{trial_num}.History;
                if strcmp(history{length(history)}, 'Abort')
                    error_abort_trials = [error_abort_trials trial_num];
                end
            end
            tl.trials([error_abort_trials]) = [];
            save('tl.mat','tl');
        end
        toc
        
        %% trace and raster plot
        % 150Hz only
        facial_mvmt = [];
        num_trials = length(tl.trials);
        time_array = [];
        time_stim_ar = [];
        
        max_trial_length = 4500;
        max_time_length = 4500;
        
        for trial_num = 1:num_trials
            
            %current trial data
            temp_time = vertcat(tl.trials{trial_num}.states.tframe)';
            temp_trial = Orofacial.Mtrace_face_split{trial_num};
            
            %check timing of each trial's stimulus onset
            temp_end = temp_time(1);
            time_stim = tl.trials{trial_num}.states(2).tframe;
            time_stim = time_stim -temp_end;
            time_stim_ar(trial_num,:) = time_stim(1);
            
            %Create a vector to capture all trials
            facial_mvmt(trial_num,:) = [temp_trial nan(1,max_trial_length-length(temp_trial))];
            time_array(trial_num,:) = [temp_time nan(1,max_time_length-length(temp_time))];
            
            % Get startle on same basis
            if ~isempty(tl.trials{trial_num}.analog.sts)
                startle(trial_num,:) = tl.trials{trial_num}.analog.data(...
                    tl.trials{trial_num}.analog.sts(2)- 10000 : ...
                    tl.trials{trial_num}.analog.sts(2) + 10000,1);
                startle_filt(trial_num,:) = filter_abr(startle(trial_num,:)', 10000, 64, 10, 2000);
            else
                startle(trial_num,:) = nan(1,20001);
                startle_filt(trial_num,:) = nan(1,20001);
            end
            
        end
        
        %% processing data to account for baseline
        % facial_mvmt(facial_mvmt>20) = NaN; % Elimintate outliers
        facial_mvmt_z = facial_mvmt;
        facial_mvmt_z(~isnan(facial_mvmt)) = zscore(facial_mvmt(isnan(facial_mvmt)~=1),[],'all');
        
        %% Plot Facial movement
        
        % get rid of the first block
        facial_mvmt_z(1:sum([tl.params.Tosca.Schedule.Families.Number]),:) = [];
        time_array(1:sum([tl.params.Tosca.Schedule.Families.Number]),:) = [];
        time_stim_ar(1:sum([tl.params.Tosca.Schedule.Families.Number]),:) = [];
        
        fme_vec = zeros(height(facial_mvmt_z), max_trial_length+1); %hind foot vec for interpolated data, +1 so that plot x,y are same length
        
        %plot trial by trial motion energy against time
        for t = 1:height(facial_mvmt_z)
            %gather the non nan trial data (time and amplitude)
            temp_x = time_array(t,:);
            temp_x = temp_x(~isnan(temp_x)); %cannot interpolate nans
            temp_x = temp_x-temp_x(1);%start at 0s
            %align times by sound onset
            sound_onset = time_stim_ar(t);
            time_diff = sound_onset - 1;
            frame_diff = round(time_diff*150);
            
            temp_y = facial_mvmt_z(t,:);
            temp_y = temp_y(~isnan(temp_y));%cannot interpolate nans
            
            %interpolate so that fs = 150 hz
            try
                int_y = interp1(temp_x,temp_y, temp_x(1):1/150:temp_x(end));
            catch errorinfo
                disp('error with avi, could not interpolate');
                int_y = temp_y;
            end
            %filter and align by stage two trial time (sound onset time)
            if nanstd(int_y(1:149)) > .95 %if there is high variability before sound comes on
                fme_vec(t,:) = nan;
            else
                if frame_diff == 0
                    fme_vec(t,:) = [int_y nan(1, (max_trial_length+1)-length(int_y))];
                elseif frame_diff < 0
                    fme_vec(t,:) = [nan(1,abs(frame_diff)) int_y nan(1, (max_trial_length+1)-(length(int_y) + abs(frame_diff)))];
                elseif frame_diff > 0
                    fme_vec(t,:) = [int_y(frame_diff+1: end) nan(1, (max_trial_length+1)-length(int_y(frame_diff+1: end)))];
                end
                fme_vec(fme_vec>20) = NaN;
            end
        end
        
        %raster plot
        % figure;
        % imagesc(fme_vec)
        % set(gca,'XTick',(0:15:4500))
        % set(gca,'XTickLabel', (0:15:4500)./150)
        % xlabel('Time re Sound Onset')
        % ylabel('Facial Movement Amplitude')
        % title(['Raster Plot-- Mouse ', num2str(mouse_num),' Sess ', num2str(session_num),' Run ', num2str(run_num)])
        % xl = xline(150,'--',{'Sound Onset'});
        % xl.LabelHorizontalAlignment = 'Left';
        
        %% separate facial movement and startle into groups
        
        %get rid of first block
        startle(1:sum([tl.params.Tosca.Schedule.Families.Number]),:) = [];
        tl.trials(1:sum([tl.params.Tosca.Schedule.Families.Number])) = [];
        
        %a little renaming and zscoring
        for t = 1:height(fme_vec)
            nbn_mat(t,:) = fme_vec(t,:);
            nbn_freq(t,:) = tl.trials{t}.State_2.Sound.Filter.CF;
            nbn_level(t,:) = tl.trials{t}.State_2.Sound.Level.dB_SPL;
            nbn_s(t,:) = startle(t,:);
        end
        nbn_startle = nbn_s;
        nbn_startle(~isnan(nbn_s)) = zscore(nbn_s(isnan(nbn_s)~=1),[],'all');
        
        %% Plots
        startle_time = (0:size(nbn_startle,2)-1)./10000;
        
        %narrowband
        nbn_freqs = unique(nbn_freq);
        nbn_levs = unique(nbn_level);
        
        figure
        cmap = hot;
        cmap = cmap(1:.5:256,:);
        cmap = cmap(round(linspace(1,350,length(nbn_levs))),:);
        
        cmap2 = jet;
        cmap2 = cmap2(1:.5:256,:);
        cmap2 = cmap2(round(linspace(1,350,length(nbn_levs))),:);
        
        nbn_mat_z = [];
        nbn_startle_z = [];
        for f = 1:length(nbn_freqs)
            subplot(2,2,f)
            for i = 1:length(nbn_levs)
                temp_y = nanmean(nbn_mat(find(nbn_freq == nbn_freqs(f) &...
                    nbn_level == nbn_levs(i)),:));
                face_time = (0:length(temp_y)-1)./150;
                plot(face_time,temp_y,'Color',cmap(i,:))
                hold on
                %if you want to plot startle:
                %         l = plot(startle_time,nanmean(nbn_startle(find(nbn_freq == nbn_freqs(f) &...
                %             nbn_level == nbn_levs(i)),:)),'--','Color',cmap2(i,:));
                %         lm(i) = l(1);
                %         hold on
                fra_n(f,i) =  nanmean(temp_y(150:195));
                nbn_mat_z = [nbn_mat_z; nbn_mat(find(nbn_freq == nbn_freqs(f) &...
                    nbn_level == nbn_levs(i)),:)];
                nbn_startle_z = [nbn_startle_z; nbn_startle(find(nbn_freq == nbn_freqs(f) &...
                    nbn_level == nbn_levs(i)),:)];
            end
            box off
            title(append(string(nbn_freqs(f)), ' kHz'))
            xlim([0 2])
            ylim([-1 10])
            set(gca,'XTick',0:.1:2)
            set(gca,'XTickLabel',-1:.1:1)
            xlabel('Time re: sound onset(s)')
            ylabel('Facial Movement (z-scored)')
            xl = xline(1,'--',{'Sound','Onset'});
            xl.LabelHorizontalAlignment = 'Left';
            legend(append(string(nbn_levs),' dB SPL'))
            %     legend(lm, append(string(nbn_levs),' dB SPL')) %startle legend
        end
        sgtitle(['Facial Movement: 1 oct Bandwidth NBN Responses'])
        
        subplot(2,2,[3:4])
        cmap = jet;
        cmap = cmap(1:60:end,:);
        for f = 1: size(fra_n,1)
            plot(nbn_levs,fra_n(f,:),'-o','Color',cmap(f,:))
            hold on
        end
        xlim([25 105])
        set(gca,'XTick',30:10:100)
        xlabel('Intensity (dB SPL)')
        ylabel('z-scored facial motion')
        box off
        lgd = legend(append(string(nbn_freqs), ' kHz'),'Location','northwest');
        title(lgd,'Center Frequency')
        title('Facial Movement vs. Intensity')
        
        %%
        cd(savingDir);
        
        %Narrowband noise
        eval(sprintf('%s_chat%02d_mat_z = nbn_mat_z;', baseline_type, mouse_num));
        eval(sprintf('%s_chat%02d_startle_z = nbn_startle_z;', baseline_type, mouse_num));
        save(sprintf('%s_chat%02d_mat_z.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_mat_z', baseline_type, mouse_num));
        save(sprintf('%s_chat%02d_startle_z.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_startle_z', baseline_type, mouse_num));
    end
end

% %Narrowband noise
% eval(sprintf('%s_chat%02d_nbn_mat = nbn_mat;', baseline_type, mouse_num));
% eval(sprintf('%s_chat%02d_nbn_freq = nbn_freq;', baseline_type, mouse_num));
% eval(sprintf('%s_chat%02d_nbn_level = nbn_level;', baseline_type, mouse_num));
% eval(sprintf('%s_chat%02d_nbn_startle = nbn_startle;', baseline_type, mouse_num));
% save(sprintf('%s_chat%02d_nbn_mat.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_nbn_mat', baseline_type, mouse_num));
% save(sprintf('%s_chat%02d_nbn_freq.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_nbn_freq', baseline_type, mouse_num));
% save(sprintf('%s_chat%02d_nbn_level.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_nbn_level', baseline_type, mouse_num));
% save(sprintf('%s_chat%02d_nbn_startle.mat', baseline_type, mouse_num), sprintf('%s_chat%02d_nbn_startle', baseline_type, mouse_num));
