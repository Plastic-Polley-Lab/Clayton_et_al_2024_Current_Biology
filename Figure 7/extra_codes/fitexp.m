function model = fitexp(lags, normalizedACF, fig)
options = fitoptions('exp1');
options.StartPoint = [1 -0.5];
% options.Upper = [Inf 0];
options.Upper = [1 0];
options.Lower = [1 -Inf];
[curve,gof] = fit(lags',normalizedACF','exp1',options);
model.tau = -1/(curve.b);
model.tau_gof_single = gof.adjrsquare;
model.exp_fit_type = 1;
if (gof.adjrsquare<0.75)
    model.exp_fit_type = 2;
    options = fitoptions('exp2');
    options.StartPoint = [1 -0.5 0.5 -0.3];
    options.Upper = [Inf 0 Inf -0.01];
    %     options.Lower = [0 -Inf 0 -Inf];
    [curve,gof] = fit(lags',normalizedACF','exp2',options);
    tau1 = -1/(curve.b);
    tau2 = -1/(curve.d);
    model.tau_final = (curve.a*tau1+curve.c*tau2)/(curve.a+curve.c);
    if model.tau_final>300
        model.tau_final = model.tau;
        model.exp_fit_type = 1;
    end
else
    model.tau_final = model.tau;
end
model.tau_gof_single = gof.adjrsquare;
model.f = curve;




if fig ==1
    figure
    h = scatter(lags, normalizedACF, [], 'ko', 'filled');
    hold on
    h1= plot(model.f);
    %         set(h1, 'Color', '-r')
    %     legend([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10)], ...
    %         {sprintf('tau = %2.2f ms', model(1).tau), sprintf('tau = %2.2f ms', model(2).tau), ...
    %         sprintf('tau = %2.2f ms', model(3).tau), sprintf('tau = %2.2f ms', model(4).tau),...
    %         sprintf('tau = %2.2f ms', model(5).tau), sprintf('tau = %2.2f ms', model(6).tau), ...
    %         sprintf('tau = %2.2f ms', model(7).tau), sprintf('tau = %2.2f ms', model(8).tau), ...
    %         sprintf('tau = %2.2f ms', model(9).tau), sprintf('tau = %2.2f ms', model(10).tau)})
    legend([h,h1], {'data', sprintf('tau = %2.2f ms', model.tau_final)})
    xlabel('Time (s)')
    ylabel('ACF')
    set(gcf,'position',[100,200,400,400])
%     close
end
% save the result
