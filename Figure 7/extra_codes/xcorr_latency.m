function [normalizedACF_p, normalizdACF_xcorr, lags] = xcorr_latency(firingRate, t, bin, fig)



 
numLags = length(firingRate) -1;
firingRate_smooth = smooth(firingRate,3);

if fig ==1
    figure
    subaxis(1, 2, 1, 'sh', 0.1, 'sv', 0.1)
    plot(t, firingRate_smooth,  '-k')
    ylabel('Change in FR (Hz)')
    xlabel('Time (ms)')
end

% using the auto-correlation
[normalizedACF_p, lags] = autocorr(firingRate_smooth, numLags);
% using the xcorr
[c,lags_temp] = xcorr(firingRate_smooth,'coeff');
lags = lags_temp((length(lags_temp)+1)/2:end);
normalizedACF = c((length(lags_temp)+1)/2:end);
normalizdACF_xcorr = normalizedACF;

if fig ==1
    subaxis(1,2,2,'sh', 0.1, 'sv', 0.1)
    h(1) = plot(lags*bin, normalizedACF, '-k')
    hold on
    h(2) = plot(lags*bin, normalizedACF_p, '-b')
    plot([0,100], [0,0], '--', 'color', [0.5, 0.5, 0.5])
    legend([h(1), h(2)], 'xcorr', 'autocorr')
    axis on
    ylim([-0.4,1])
    ylabel('ACF')
    xlabel('Time (ms)')
end
% 
% 
% % for i = 1:2  
% %     idx = find(normalizedACF(i,:)< 0);
% %     lags_fit{i} = lags(1:idx(1));
% %     normalizedACF_fit{i} = normalizedACF(i,1:idx(1));
% %     % title(['Jitter ascending; Unit', num2str(neuron_num)])
% %     model(i).f = fit(lags_fit{i}',normalizedACF_fit{i}','exp1');
% %     model(i).tau = -1/model(i).f.b;
% % end
% 
% 
%     options = fitoptions('exp1');
%     options.StartPoint = [1 -0.5];
%     % options.Upper = [Inf 0];
%     options.Upper = [1 0];
%     options.Lower = [1 -Inf];
%     [curve,gof] = fit(lags(i,:)',normalizedACF(i,:)','exp1',options);
%     model(i).tau = -1/(curve.b);
%     model(i).tau_gof_single = gof.adjrsquare;
%     model(i).exp_fit_type = 1;
%     if (gof.adjrsquare<0.75)
%         model(i).exp_fit_type = 2;
%         options = fitoptions('exp2');
%         options.StartPoint = [1 -0.5 0.5 -0.3];
%         options.Upper = [Inf 0 Inf -0.01];
%         %     options.Lower = [0 -Inf 0 -Inf];
%         [curve,gof] = fit(lags(i,:)',normalizedACF(i,:)','exp2',options);
%         tau1 = -1/(curve.b);
%         tau2 = -1/(curve.d);
%         model(i).tau_final = (curve.a*tau1+curve.c*tau2)/(curve.a+curve.c);
%         if model(i).tau_final>100
%             model(i).tau_final = model(i).tau;
%             model(i).exp_fit_type = 1;
%         end
%     else
%         model(i).tau_final = model(i).tau;
%     end
%     model(i).tau_gof_single = gof.adjrsquare; 
%     model(i).f = curve;





% if fig ==1
%     subaxis(2,2,4,'sh', 0.1, 'sv', 0.1)
%     hold off
%     for i = 1:2
%         h(i) = scatter(lags(i,:), normalizedACF(i,:), [], CT(i,:), 'filled');
%         hold on
%         h1= plot(model(i).f);
%         set(h1, 'Color', CT(i,:))
%     end
% %     legend([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10)], ...
% %         {sprintf('tau = %2.2f ms', model(1).tau), sprintf('tau = %2.2f ms', model(2).tau), ...
% %         sprintf('tau = %2.2f ms', model(3).tau), sprintf('tau = %2.2f ms', model(4).tau),...
% %         sprintf('tau = %2.2f ms', model(5).tau), sprintf('tau = %2.2f ms', model(6).tau), ...
% %         sprintf('tau = %2.2f ms', model(7).tau), sprintf('tau = %2.2f ms', model(8).tau), ...
% %         sprintf('tau = %2.2f ms', model(9).tau), sprintf('tau = %2.2f ms', model(10).tau)})
%     legend([h(1), h(2)], ...
%         {sprintf('reg tau = %2.2f ms', model(1).tau_final), sprintf('rand tau = %2.2f ms', model(2).tau_final)})   
%     xlabel('Time (s)')
%     ylabel('ACF')
%     set(gcf,'position',[100,200,800,800])
%     suptitle(['Neuron ', num2str(neuron_num)])
%     
%     if neuron_num < 10
%         print(['Jitter-ascending_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
% %            set(gcf, 'Color', 'w')
% %            export_fig(['Jitter-ascending_response_0', num2str(neuron_num)],  '-png')
%     else
%         print(['Jitter-ascending_response_', num2str(neuron_num)],'-dpdf','-bestfit')
% %            set(gcf, 'Color', 'w')
% %            export_fig(['Jitter-ascending_response_', num2str(neuron_num)],  '-png')
%     end
%     close
% end
% % save the result
% pattern_result.psth_summary = psth_summary;
% pattern_result.firingRate   = firingRate;
% pattern_result.normalizedACF= normalizedACF;
% pattern_result.lags = lags;
% % pattern_result.lags_fit = lags_fit;
% % pattern_result.normalizedACF_fit = normalizedACF_fit;
% pattern_result.model = model;
% end