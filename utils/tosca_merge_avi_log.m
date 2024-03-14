function TL = tosca_merge_avi_log(TL, AVI)

for k = 1:length(TL.trials)
    
    if abs(TL.trials{k}.start - AVI(k).toscaTime) > 1
        [TL, AVI] = tosca_repair_video(TL, AVI);
        %       error('Possible video error. Please tell Ken.');
    end
    
    TL.trials{k}.aviFile = AVI(k).aviFile;
    TL.trials{k}.ttosca = AVI(k).toscaTime;
    TL.trials{k}.frameRate = AVI(k).frameRate;
    
    states = TL.trials{k}.states;
    
    %Note: mismatch between number of states for AVI vs. TL
    if length(AVI(k).states) > length(states) 
        AVI(k).states = AVI(k).states(1:length(states));
        
    elseif length(AVI(k).states) < length(states)
        states = states(1:length(AVI(k).states));
    end
    
    TL.trials{k} = rmfield(TL.trials{k}, 'states');
    
    for ks = 1:length(states)
        st = states(ks);
        
        st.frames = AVI(k).states(ks).frameNum;
        st.tframe = AVI(k).states(ks).tframe;
        
        TL.trials{k}.states(ks) = st;
    end
end

% function TL = tosca_merge_avi_log(TL, AVI)
%
% for k = 1:length(TL.trials)
%
%    if abs(TL.trials{k}.start - AVI(k).toscaTime) > 1
%       [TL, AVI] = tosca_repair_video(TL, AVI);
% %       error('Possible video error. Please tell Ken.');
%    end
%
%    TL.trials{k}.aviFile = AVI(k).aviFile;
%    TL.trials{k}.ttosca = AVI(k).toscaTime;
%    TL.trials{k}.frameRate = AVI(k).frameRate;
%
%    states = TL.trials{k}.states;
%    %Note: mismatch between number of states for AVI vs. TL
%    if length(AVI(k).states) ~= length(states)
%        states = states(1:length(AVI(k).states));
%    end
%
%    TL.trials{k} = rmfield(TL.trials{k}, 'states');
%
%    for ks = 1:length(states)
%       st = states(ks);
%
%       st.frames = AVI(k).states(ks).frameNum;
%       st.tframe = AVI(k).states(ks).tframe;
%
%       TL.trials{k}.states(ks) = st;
%    end
% end


%         if(length(AVI(k).states) > length(states)) && length(states) == 3
%             AVI(k).states = AVI(k).states(1:length(states));
%
%         elseif (length(AVI(k).states) > length(states)) && length(AVI(k).states) == 3
%             states = states(1:length(AVI(k).states));
%
%         elseif (length(AVI(k).states) < length(states)) && length(states) == 3
%             AVI(k).states = AVI(k).states(1:length(states));
%
%         elseif (length(AVI(k).states) < length(states)) && length(AVI(k).states) == 3
%             states = states(1:length(AVI(k).states));
%         else
%         end