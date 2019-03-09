y1 = @(x) exp(x) / 4.5;
fplot(y1, [0.5, 2.5])
hold on
y2 = @(x) (exp(1.5) * (x - 1.5) + exp(1.5))/4.5;
fplot(y2, [0.6, 2.5])
y3 = @(x) y2(x) - 0.5
fplot(y3, [0.7, 2.5])
axis([0, 3, -0.5, 3])
scatter([1.5, 2.3, 2.3, 1, 1], [exp(1.5) / 4.5, y1(2.3), y2(2.3), y1(1), y2(1)], '.k')
plot([1, 1], [-0.5, y1(1)], '--')
plot([0.95, 1], [y1(1), y1(1)], '-k')
plot([0.95, 1], [y2(1), y2(1)], '-k')
plot([0.95, 1], [y3(1), y3(1)], '-k')
plot([2.3, 2.3], [-0.5, y1(2.3)], '--')
plot([2.25, 2.3], [y1(2.3), y1(2.3)], '-k')
plot([2.25, 2.3], [y2(2.3), y2(2.3)], '-k')
plot([2.25, 2.3], [y3(2.3), y3(2.3)], '-k')

annotation('doublearrow', [2.25 2.25] / 3.135, [y1(2.3) + 0.01, y2(2.3) + 0.145] / 3, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [2.25 2.25] / 3.135, [y2(2.3) + 0.145, y3(2.3) + 0.29] / 3, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [0.95 0.95] / 2.485, [y1(1) + 0.13, y2(1) + 0.185] / 2, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [0.95 0.95] / 2.485, [y2(1) + 0.185, y3(1) + 0.455] / 2, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')

text(2.5994, 1.9871, '$l$', 'FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(2.5994, 1.9871 - 0.5, '$l_{1}$', 'FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.6, 0.9638, '$V_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

text(2.2, 2.3376, '$V_{A1}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(2.35, 1.6998, '$V_{B1}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(0.95, 0.7325, '$V_{A2}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.05, 0.4030, '$V_{B2}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

text(2.2, (y1(2.3) + y2(2.3))/2, '$\delta$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(2.2, (y3(2.3) + y2(2.3))/2, '$t_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(0.9, (y1(1) + y2(1))/2, '$\delta$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(0.9, (y3(1) + y2(1))/2, '$t_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$V$','FontName','Times New Roman', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'YTick', []);