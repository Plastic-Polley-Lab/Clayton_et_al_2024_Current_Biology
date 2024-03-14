%%
function cp = extract_cp(Onset, Offset, H, sign)
% sign = 1;
% Onset  = wt_psth_pop.exc_Time_On;
% Offset = wt_psth_pop.exc_Time_Off;
% H      = wt_psth_pop.H;
for i = 1:length(Onset)
    if isnan(Onset{i})
        cp.firstOnset(i)  = NaN;
        cp.firstH(i)      = NaN;
        cp.firstOffset(i) = NaN;
        cp.lastOffset(i)  = NaN;
        cp.duration(i)    = NaN;
    else
        temp =      Onset{i};
        temp_offset = Offset{i};
        temp_H      = H{i};
        indx = find(temp_H == sign);
        if isempty(indx)
            warning('something is absurd')
            fprintf('Checking Neuron # %d \n', i)
            cp.firstOnset(i)  = NaN;
            cp.firstH(i)      = NaN;
            cp.firstOffset(i) = NaN;
            cp.lastOffset(i)  = NaN;
            cp.duration(i)    = NaN;
        else
            iter = [];
            for j = indx(1):length(temp_H)
                if temp_H(j) == sign
                    iter = [iter,j];
                else
                    break
                end
            end
            z_indx = [];
            if length(iter)==1
               z_indx = 1;
            else
                for z = 1:(length(iter)-1)
                    if temp(iter(z+1)) == temp_offset(iter(z)) % find the continuous responses
                       z_indx = z + 1;
                    else
                       z_indx = z;
                        break
                    end
                end
            end
            
            cp.firstOnset(i)  = temp(iter(1));
            cp.firstOffset(i) = temp_offset(iter(1));
            cp.firstH(i)      = temp_H(iter(1));
            cp.lastOffset(i)  = temp_offset(iter(z_indx));
            cp.duration(i)    = cp.lastOffset(i) - cp.firstOnset(i);
        end
    end
end