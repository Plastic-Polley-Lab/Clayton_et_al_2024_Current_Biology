function [frqTuning lvlTuning cFRA spontDistr fraDistr dprime FRAmask] = CleanFRA_KenGuo(FRA, mask, param)

%FRA should be x by x 
%FRA = FRA'; % NOTE Needs level x frequency as input to work 
%Updated this 10.17.23

% Pad the FRA to eliminate edge effects from median filter
padFRA = zeros(size(FRA,1)+2,size(FRA,2)+2);
padFRA(2:end-1,2:end-1) = FRA;
padFRA(:,1) = padFRA(:,2);
padFRA(:,end) = padFRA(:,end-1);
padFRA(1,:) = padFRA(2,:);
padFRA(end,:) = padFRA(end-1,:);

% Median filter (x4)
padFRA = medfilt2(padFRA,[4 4]);
padFRA = medfilt2(padFRA,[3 3]);
padFRA = medfilt2(padFRA,[3 3]);
padFRA = medfilt2(padFRA,[3 3]);

% Extract central portion of median filtered FRA
mfFRA = padFRA(2:end-1,2:end-1);

% Compute frequency tuning function
frqTuning = sum(mfFRA); % sum across levels
frqTuning = medfilt1(frqTuning, 2); % N-point moving median
% Pad to handle edge effects from filtfilt
% -- why 0.5?
frqTuning = [ones(1,5)*frqTuning(1)*0.5 frqTuning ones(1,5)*frqTuning(end)*0.5];
frqTuning = filtfilt(ones(1,3)/3, 1, frqTuning); % 3-point boxcar filter, forward and backward

% Compute level tuning function
lvlTuning = sum(mfFRA,2)'; % sum across frequency note there was an errant transpose here. 
% Pad to handle edge effects from median filter
% -- why asymmetric?
lvlTuning = [ ...
   ones(1,5)*mean(lvlTuning(1:2))*0.5 ...
   lvlTuning ...
   ones(1,2)*lvlTuning(end) ...
   ];
% Median filter (x2)
lvlTuning = medfilt1(lvlTuning,4);
lvlTuning = medfilt1(lvlTuning,4);


% Level threshold criterion is mean rate + some multiple of the s.d. of the
% first 6 points (only the last of which is an actual data point?!?)
lTThr = mean(lvlTuning(1:6)) + std(lvlTuning(1:6))*(param*10);
% Find last point where rate is below threshold criterion
idx = find(lvlTuning<=lTThr, 1, 'last');
% Subtract 5: to account for extraction of central portion
lTThrInd = max(idx-5, 1);
lvlTuning = lvlTuning(6:end-2); % extract central portion of median-filtered data

% 
sFT = sort(frqTuning);
% Frequency threshold criterion is mean + s.d. of the 6th-9th largest
% rates?!?
fTThr = mean(sFT(6:9)) + std(sFT(6:9));
[~, peakInd] = max(frqTuning); %

% Low frequency side: find last point before the peak below rate criterion
for i = peakInd:-1:1;
   if (frqTuning(i) <= fTThr)
      break;
   end;
end;
startInd = max(1, i-5); % account for filter padding

% High frequency side: find first point beyond peak below rate criterion
for i = peakInd:length(frqTuning);
   if (frqTuning(i) <= fTThr)
      break;
   end;
end;
endInd = min(size(FRA,2), i); % account for filter padding

frqTuning = frqTuning(5:end-4); % extract central portion of filtered data

if isempty(mask),
    FRAmask = zeros(size(FRA)); % initialize FRA mask
    if endInd - startInd > 1, %This sets the threshold for frequency tuning, must be at least 1 octave 
        fTmask = zeros(size(frqTuning));
        fTmask(startInd:endInd) = frqTuning(startInd:endInd) - fTThr; %Set values to 0 point 
        fTmask = length(lvlTuning) - round((fTmask/max(fTmask))*(length(lvlTuning)-lTThrInd + 1)) + 1;
        if sum(fTmask<=0)>0
            disp('invalid mask')
        else
            for i = 1:size(FRA,2);
                FRAmask(fTmask(i):end,i) = 1; %For values above threshold at each frequency, assign threshold
            end
        end
    end
else
    FRAmask = mask'; %Retranspose
end

% figure(1)
% imagesc(FRAmask) 

cFRA = FRAmask.*FRA;


fraList = FRA(FRAmask > 0);
spontList = FRA(FRAmask < 1);
FRAmaskupper = FRAmask;
FRAmaskupper(1:round(size(FRAmask,1)/2),:) = 1;
spontList2 = FRA(FRAmaskupper < 1);
randTrialNum = 1000;
sampNum = 40;
spontDistr = zeros(randTrialNum,1);
fraDistr = zeros(randTrialNum,1);
for i = 1:randTrialNum;
   if (length(fraList)*0.5 > sampNum)
      fraDistr(i) = mean(fraList(ceil(rand(sampNum,1)*(length(fraList)-1))+1));
      spontDistr(i) = mean(spontList(ceil(rand(sampNum,1)*(length(spontList)-1))+1));
   else
      %         disp(spontList2);
      fraDistr(i) = mean([fraList(ceil(rand(round(length(fraList)*0.5),1)*(length(fraList)-1))+1);...
         spontList2(ceil(rand(sampNum-round(length(fraList)*0.5),1)*...
         (length(spontList2)-1))+1)]);
      spontDistr(i) = mean(spontList(ceil(rand(sampNum,1)*(length(spontList)-1))+1));
   end;
end;

dprime = (mean(fraDistr)-mean(spontDistr))*2/(std(fraDistr)+std(spontDistr));
% dprime = (mean(fraDistr)-mean(spontDistr))/(std(fraDistr));


frqTuning = sum(cFRA);
lvlTuning = sum(cFRA,2)';
