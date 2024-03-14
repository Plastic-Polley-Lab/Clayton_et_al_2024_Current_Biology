% single neuron decoding based on psth_template
% here I'll use leave-one-out as the way of doing cross-validation
% try to split the data into training set and testing set
%%
function options = single_neuron_decoding_PSTH_v2(data, options)
% data = summary(12); % data for decoding
trig=fieldnames(data);
nev=numel(trig);
%% Set the parameters for decodiing: windows and binsize; 
windows    = [251, 1000]; % post-stimulus window: stimuli starts at 251 ms after 
binsize    = 1 ; % 1 ms
smooth_set = 1;  % whether smooth the PSTH before decoding
%% re-raster the psth if needed
if options. binsize == 1
    psth = data;
    options.windows = windows;
    options.binsize = binsize;
else
    for i = 1:nev % re-raster the spike with different binsize
        for j = 1:length(data.(trig{i}).raster)
            spike_times = data.(trig{i}).raster(j).ts;
            spikes = histcounts(spike_times,0:binsize:windows(2));
            if isempty(spikes) % no spikes in that trial
                raster(j,:) = zeros(length(0:binsize:windows(2))-1,1); % pad zeros
            else
                raster(j,:) = spikes;
            end
        end
        psth.(trig{i}).scmatrix = raster;
        clear raster
    end
    t = (0:binsize:windows(2)) + binsize/2 ; % get the time of the psth
    t = t(1:end-1);
    options.windows(1) = max(find(t<windows(1))); % re-get the windows after re-raster
    options.windows(2) = length(t);
    options.binsize    = binsize;
end
%% get the total numebr of classes

% assuming you have the same number of trials/or
% find smallest number of trials across events
ntrials=zeros(1,numel(trig));
for E=1:numel(trig)
     ntrials(E)=size(psth.(trig{E}).raster,2);
end
%%
options.mintrials=min(ntrials); % pick smallest number of trials across conditions for decoding\
options.nclass = nev;
options.nrun   = nev * min(ntrials);
options.trueClass = zeros(options.nrun,1);
options.Classified = zeros(options.nrun,1);



for i = 1:nev
    testing{i} = randperm(options.mintrials);
end
clear temp
for i = 1:nev
    for j = 1:length(testing{i})
        temp(:, j) = setdiff(1:options.mintrials, testing{i}(j))'; % get the training dataset for each testing
    end
    training{i} = temp;
    clear temp
end

% start the decoding
iter = 1;
for i = 1: length(testing{1}) % test each test set
    for j = 1: options.nclass
        temp.(trig{j}) = mean(psth.(trig{j}).scmatrix(training{j}(:,i), options.windows(1):options.windows(2)),1);
        test.(trig{j}) = psth.(trig{j}).scmatrix(testing{j}(i), options.windows(1):options.windows(2));
        if smooth_set
            temp.(trig{j}) = smoothts(temp.(trig{j}), 'g', 5, 1); % smooth the traces
            test.(trig{j}) = smoothts(test.(trig{j}), 'g', 5, 1); % smooth the traces
        else
        end
    end
    % calculate the euclidean distance of each test to all other class
    for j = 1:options.nclass
        for z = 1:options.nclass
            distance.(trig{j})(z) = sqrt(sum((test.(trig{j}) - temp.(trig{z})).^2));
        end
        [M, I] = min(distance.(trig{j}));
        options.trueClass(iter) = j;
        options.Classified(iter) = I;
        iter = 1 + iter;
    end
    
end

options.confMatrix = rot90(confusionmat(options.trueClass,options.Classified));
% figure;
% imagesc(options.confMatrix./options.mintrials);
% xlabel('True Classes')
% ylabel('Classified Classes')
% colormap(hot)

for i = 1:nev % get the decoding performance for each class
    idx = find(options.trueClass == i);
    idx_true = find(options.Classified(idx) == i);
    options.performance(i) = length(idx_true)/length(idx);
end
options.avg_performance = length(find(options.trueClass == options.Classified))/length(options.Classified);
