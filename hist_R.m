hist(R, 40)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 105 185] / 255;
h.EdgeColor = 'k';
title('指数增长模型', 'FontSize', 20)
%title('二次多项式模型', 'FontSize', 20)