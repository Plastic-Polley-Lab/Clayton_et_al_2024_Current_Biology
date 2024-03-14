% add by ke to slight disrupt the regular rhythm
function t_jitter = interval_jitter(intervals, reg_interval, jitters, Variable, jitter_type)
if jitters ~=0
    CycleLength = Variable.CycleLength;
    lowB = 100;  % minimal interval
    highB = 1000; % add be Ke, increase this value very high to remove this constrain. 
    t_jitter = []; % within each frame; the interval goes up/down  a fixed value within certain range; but keep the frame duration same
    switch jitter_type
        case 'jitter-fixed'
            for k = 1:Variable.NumCycles
                jitter_time = randi(floor(mean(intervals) * jitters));
                seq = repmat(jitter_time, size(reg_interval));
                seq(1:CycleLength/2) = -jitter_time;
                outboundary = 1;
                while outboundary
                    temp = seq(randperm(CycleLength));
                    temp = reg_interval + temp;
                    min_idx = find(temp < lowB);
                    max_idx = find(temp > highB);
                    if isempty(min_idx) && isempty(max_idx)
                        outboundary =0;
                    end
                end
                t_jitter = [t_jitter, temp];
            end
    
        case 'jitter-rand'
            t_jitter = []; % within each frame; change the interval by a random value within certain range; but keep the frame duration same
            for k = 1:Variable.NumCycles
                temp = reg_interval;
                jitter_time = floor(mean(intervals) * jitters);
                outboundary = 1;
                while outboundary
                    jitters_interval = randomSeq(-jitter_time:jitter_time, CycleLength);
                    temp = reg_interval + jitters_interval;
                    min_idx = find(temp < lowB);
                    max_idx = find(temp > highB);
                    if isempty(min_idx) && isempty(max_idx)
                        outboundary =0;
                    end
                end
                t_jitter = [t_jitter, temp];
            end
            
        case 'jitter-scale'
            t_jitter   = []; % within each frame; keep the ratio between interval the same; but scale up and down the frame duration;
            jitter_seq = -jitters : 0.01 : jitters; % modified by Ke to make the jitter bigger
%             jitter_seq = [-jitters: 0.01: -(jitters - 0.1), (jitters - 0.1):0.01:jitters];
            temp = reg_interval;
            for k = 1:Variable.NumCycles
                rng('shuffle')
                jitter_scale = jitter_seq(randi(length(jitter_seq)));
                t_jitter = [t_jitter, temp + temp*jitter_scale];
            end
%         case 'jitter-swap'
%             swap_idx = randperm(length(reg_interval),Variable.swap); % swap interval
%             swap_interval = reg_interval;
%             swap_interval(swap_idx) = reg_interval(flip(swap_idx));
%             t_swap = repmat(swap_interval, 1, Variable.NumCycles);
    end

    
else
    error('This program is for jitting the intervals')
end
