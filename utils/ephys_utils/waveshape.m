

function [pt_delay_ms,unit_temp,norm_temp]=waveshape(stTemplate,templateMin)
%stTemplate = waveforms across channels
%templateMin = unit location 

        unit_temp = stTemplate(templateMin,:);

        %normalize the unit template
        norm_temp =unit_temp./max(abs(unit_temp));

        %grab statistics about peak/trough
        [tVal,tIdx] = min(unit_temp);
        [pVal,pIdx] = max(unit_temp(tIdx:end));
        pIdx = pIdx+(tIdx-1);
        pAmp = pVal;
        pt_delay = pIdx - tIdx;
        pt_ratio = pVal/abs(tVal); %this is what jennifer did
        Fs = 24414;
        pt_delay_ms = (pt_delay./Fs).*1000;
 
end 
