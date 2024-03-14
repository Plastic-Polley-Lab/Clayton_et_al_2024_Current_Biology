function [cFRA spontDistr fraDistr dprime FRAmask] = dprime_FRA(FRA)
%     FRA =temp(2).fra;
    mfFRA = mean(FRA,3);

    frqTuning = sum(mfFRA,2);
%     frqTuning = medfilt1(frqTuning,5);
%     frqTuning = [ones(1,5)*frqTuning(1)*0.5,frqTuning,ones(1,5)*frqTuning(end)*0.5];
%     frqTuning = filtfilt(ones(1,3)/3, 1, frqTuning);

    lvlTuning = mean(mfFRA,1);
    % get the baseline
    lvl_lowest = mean(squeeze(FRA(:,1,:)),2);
    lvl_lowest_m = mean(lvl_lowest(:));
    lvl_lowest_std = std(lvl_lowest(:));
    
    % param would be from 0 to 1
    param = 0.1;
    lTThr = lvl_lowest_m + lvl_lowest_std * (param*10);
%     lvlTuning = sum(mfFRA,2)';
%     lvlTuning = [ones(1,5)*mean(lvlTuning(1:2))*0.5,...
%         lvlTuning,ones(1,2)*lvlTuning(end)];
%     lvlTuning = medfilt1(lvlTuning,4);
%     lvlTuning = medfilt1(lvlTuning,4);

    

    sLT = sort(lvlTuning);
    lTThr = mean(sLT(1:3)) + std(sLT(1:3))*(param*10); 
    for i = length(lvlTuning):-1:1;
        if (sLT(i) <= lTThr)
            break;
        end;
    end;
    lTThrInd = max(i,1);
%     lvlTuning = lvlTuning(6:end-2);

    sFT = sort(frqTuning);
    fTThr = mean(sFT(1:3)) + std(sFT(1:3))*1;
    [tilda,peakInd] = max(frqTuning);
    for i = peakInd:-1:1;
        if (frqTuning(i) <= fTThr)
            break;
        end;
    end;
    startInd = max(1,i);
    for i = peakInd:length(frqTuning);
        if (frqTuning(i) <= fTThr)
            break;
        end;
    end;
    endInd = min(size(mfFRA,1),i);
%     frqTuning = frqTuning(5:end-4);
mask =[];
if (isempty(mask))
    FRAmask = zeros(size(mfFRA));
    if (endInd-startInd) > 1
        fTmask = zeros(size(frqTuning));
        fTmask(startInd:endInd) = frqTuning(startInd:endInd)-fTThr;
%         if ~isempty(find(fTmask<0))
%             error('check here', 'something is not working')
%         end
        fTmask = length(lvlTuning) - round((fTmask/max(fTmask))*(length(lvlTuning)-lTThrInd+1))+1;
        for i = 1:size(mfFRA,1);
            FRAmask(i, fTmask(i):end) = 1;
        end;
    end;
else
    FRAmask = mask;
end


cFRA = FRAmask.*mfFRA;

fraList = mfFRA(FRAmask > 0);
spontList = mfFRA(FRAmask < 1);
% FRAmaskupper = FRAmask;
% FRAmaskupper(1:round(size(FRAmask,1)/2),:) = 1;
% spontList2 = mfFRA(FRAmaskupper < 1);
randTrialNum = 1000;
sampNum = 10;
spontDistr = zeros(randTrialNum,1);
fraDistr = zeros(randTrialNum,1);
if ~isempty(fraList)
    for i = 1:randTrialNum
        fraDistr(i) = mean(datasample(fraList, sampNum, 'Replace', true));
        spontDistr(i) = mean(datasample(spontList, sampNum, 'Replace', true));
        
    end
    dprime = (mean(fraDistr)-mean(spontDistr))/sqrt(1/2*((std(fraDistr))^2+(std(spontDistr))^2));
else
    dprime = 0;
end
% dprime = (mean(fraDistr)-mean(spontDistr))/(std(fraDistr));

% frqTuning = sum(cFRA);
% lvlTuning = sum(cFRA,2)';
