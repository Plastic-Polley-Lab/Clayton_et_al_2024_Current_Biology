function Cost = rlf_modelCost(x, level, FRdata)

FR_model = (x(1) * exp(-(level - x(2)).^2/(2 * x(3).^2)) + x(4)).*(level<=x(2)) +  ( x(1) * exp(-(level - x(2)).^2/(2 * x(5).^2)) + x(6)).*(level>x(2));

Cost = sum((FRdata - FR_model).^2);
