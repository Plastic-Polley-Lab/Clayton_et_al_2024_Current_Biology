function [distance_ap, pupil_size] = DLC_pupil(data)

% get the distance based on anterior and posterior marker of pupil
% anterior: is the max x location and posterior is the min x location
for i = 1:length(data)
    pupil_data = data(i); % loop through each frame
    [~, ind_anterior] = max(pupil_data.x);
    [~, ind_posterior] = min(pupil_data.x);
    % get the pupil 'diameter' using the anterior, posterior markers
    distance_ap(i) = sqrt((pupil_data.x(ind_anterior) - pupil_data.x(ind_posterior)).^2 + ...
        (pupil_data.y(ind_anterior) - pupil_data.y(ind_posterior)).^2);
    
    cutoff = 0.7; % frame with marker likelihood bigger than 0.7 was used, others will be interpolated
    indx = find(pupil_data.likelihood<0.7);
    if length(indx)>=3 % at lease frame less than 3 marker are less confident
        fprintf('Frame # %d has likelihood less than %1.1f \n', i, cutoff)
        pupil_size.pupil_area(i) = NaN;
        pupil_size.pupil_param{i} = NaN;
        pupil_size.total_error(i) = NaN;
    else
        pupil_data.x(indx) =[];
        pupil_data.y(indx) =[];
        pupil_data.likelihood(indx) =[];
        try
            
%             figure;
%             hold on
%             scatter(pupil_data.x,pupil_data.y,'o')
            
            
            %% fit ellipse
            x = pupil_data.x ;
            y = pupil_data.y ;
            c = [mean(x), mean(y)];
%             fun = @(a) (x-c(1)).^2/a(1).^2 + (y-c(2)).^2/a(2).^2 -1;
            fun = @(a) (x-a(1)).^2/a(3).^2 + (y-a(2)).^2/a(4).^2 -1;

            x_origin = c(1);
            y_origin = c(2);
%             a0 = [10,10];
            a0 = [c, 10,10];

            options = optimset('Display','off');
            af = lsqnonlin(fun, a0, [], [], options);
            
            % hold on
%             t = linspace(0,2*pi) ;
%             x2 = c(1) + af(1)*cos(t);
%             y2 = c(2) + af(2)*sin(t);
%             x2 = af(1) + af(3)*cos(t);
%             y2 = af(2) + af(4)*sin(t);

%             plot(x2,y2, 'k--')
%             cost_Value = (x-c(1)).^2/af(1).^2 + (y-c(2)).^2/af(2).^2 -1;
            cost_Value = (x-af(1)).^2/af(3).^2 + (y-af(2)).^2/af(4).^2 -1;

%             pupil_size.pupil_area(i) = pi * af(1) *af(2);
%             pupil_size.pupil_param{i}= [c, af(1), af(2)];
            pupil_size.pupil_area(i) = pi * af(3) *af(4);
            pupil_size.pupil_param{i}= af;
            %% calculate the fitting errors
            % anterior error; fix the y position, see the x position
            %         for z = 1:length(pupil_data.x)
            %         syms a
            %         S = vpasolve((a-c(1)).^2/af(1).^2 + (data.pupil_anterior.y(i)-c(2)).^2/af(2).^2 ==1, [a], [x_origin, inf]);
            %         error(1) = abs(double(S)-data.pupil_anterior.x(i));
            % %
            %         % posterior error; fix the y position, see the x position
            %         syms a
            %         S = vpasolve((a-c(1)).^2/af(1).^2 + (data.pupil_posterior.y(i)-c(2)).^2/af(2).^2 ==1, a, [-inf,x_origin]);
            %         error(2) = abs(double(S)-data.pupil_posterior.x(i));
            %
            %         % dorsal error; fix the x position, see the y position
            %         syms a
            %         S = vpasolve((data.pupil_dorsal.x(i)-c(1)).^2/af(1).^2 + (a-c(2)).^2/af(2).^2 ==1, a, [-inf,y_origin]);
            %         error(3) = abs(double(S)-data.pupil_dorsal.y(i));
            %
            %         % ventral error; fix the x position, see the y position
            %         syms a
            %         S = vpasolve((data.pupil_ventral.x(i)-c(1)).^2/af(1).^2 + (a-c(2)).^2/af(2).^2 ==1, a, [y_origin,inf,]);
            %         error(4) = abs(double(S)-data.pupil_ventral.y(i));
            %         pupil.total_error(i) = sum(error);
        catch
            warning('There is something wrong')
            pupil.pupil_area(i) = NaN;
            pupil.pupil_param{i} = NaN;
            pupil.total_error(i) = NaN;
            continue
        end
    end
    
end
