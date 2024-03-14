function h = lineplot_error(t, data, color)
data_m = mean(data,1);
data_sem = std(data)./sqrt(size(data,1));
h = boundedline(t, data_m, data_sem, color, 'alpha')
% set(h, 'Color', [0,1,0])