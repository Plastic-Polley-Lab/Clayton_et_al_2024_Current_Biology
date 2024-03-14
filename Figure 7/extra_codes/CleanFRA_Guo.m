function [frqTuning lvlTuning cFRA spontDistr fraDistr dprime FRAmask] = CleanFRA_Guo(FRA,mask,param)
    FRA = rot90(mean(fra, 3));
    padFRA = zeros(size(FRA,1)+2,size(FRA,2)+2);
    padFRA(2:end-1,2:end-1) = FRA;
%     padFRA(:,1) = padFRA(:,2);
%     padFRA(:,end) = padFRA(:,end-1);
%     padFRA(1,:) = padFRA(2,:);
%     padFRA(end,:) = padFRA(end-1,:);
% 
%     padFRA = medfilt2(padFRA,[4 4]);
%     padFRA = medfilt2(padFRA,[3 3]);
%     padFRA = medfilt2(padFRA,[3 3]);
%     padFRA = medfilt2(padFRA,[3 3]);

    mfFRA = padFRA(2:end-1,2:end-1);

    frqTuning = sum(mfFRA);
%     frqTuning = medfilt1(frqTuning,5);
%     frqTuning = [ones(1,5)*frqTuning(1)*0.5,frqTuning,ones(1,5)*frqTuning(end)*0.5];
%     frqTuning = filtfilt(ones(1,3)/3, 1, frqTuning);

lvlTuning = sum(mfFRA,2)';
%     lvlTuning = sum(mfFRA,2)';
%     lvlTuning = [ones(1,5)*mean(lvlTuning(1:2))*0.5,...
%         lvlTuning,ones(1,2)*lvlTuning(end)];
%     lvlTuning = medfilt1(lvlTuning,4);
%     lvlTuning = medfilt1(lvlTuning,4);

    
    % param would be from 0 to 1
    param = 0.2
    sLT = sort(lvlTuning);
    lTThr = mean(sLT(1:3)) + std(sLT(1:3))*(param*10); % before edits it is from 1:6
    for i = length(lvlTuning):-1:1;
        if (sLT(i) <= lTThr)
            break;
        end;
    end;
    lTThrInd = max(i-5,1);
    lvlTuning = lvlTuning(6:end-2);

    sFT = sort(frqTuning);
    fTThr = mean(sFT(6:9)) + std(sFT(6:9))*1;
    [tilda,peakInd] = max(frqTuning);
    for i = peakInd:-1:1;
        if (frqTuning(i) <= fTThr)
            break;
        end;
    end;
    startInd = max(1,i-5);
    for i = peakInd:length(frqTuning);
        if (frqTuning(i) <= fTThr)
            break;
        end;
    end;
    endInd = min(size(FRA,2),i);
    frqTuning = frqTuning(5:end-4);

if (isempty(mask))
    FRAmask = zeros(size(FRA));
    if ((endInd-startInd) > 4)
        fTmask = zeros(size(frqTuning));
        fTmask(startInd:endInd) = frqTuning(startInd:endInd)-fTThr;
        fTmask = length(lvlTuning) - round((fTmask/max(fTmask))*(length(lvlTuning)-lTThrInd + 1)) + 1;
        for i = 1:size(FRA,2);
            FRAmask(fTmask(i):end,i) = 1;
        end;
    end;
else
    FRAmask = mask;
end


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
