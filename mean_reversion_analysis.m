% m = median(R);
% bigger_than_median = find(R >= m);
% bigger_than_median = bigger_than_median(1 : end - 1);
% bigger_follower = bigger_than_median - 1;
% bigger_diff = R(bigger_follower) - R(bigger_than_median);
% hist(bigger_diff, 20)
% h = findobj(gca,'Type','patch');
% h.FaceColor = [0 105 185] / 255;
% h.EdgeColor = 'k';
% mean(bigger_diff)
% numel(bigger_diff(bigger_diff < 0))/numel(bigger_diff)

m = median(R);
smaller_than_median = find(R <= m);
smaller_than_median = smaller_than_median(2 : end);
smaller_follower = smaller_than_median - 1;
smaller_diff = R(smaller_follower) - R(smaller_than_median);
hist(smaller_diff, 20)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 105 185] / 255;
h.EdgeColor = 'k';
mean(smaller_diff)
numel(smaller_diff(smaller_diff > 0))/numel(smaller_diff)

title('指数增长模型', 'FontSize', 20)
%title('二次多项式模型', 'FontSize', 20)