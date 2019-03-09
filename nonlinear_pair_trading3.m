y1 = @(x) log(x) * 4;
fplot(y1, [1.5, 9.5])
axis([0, 11, 0, 11])
hold on
y2 = @(x) (1/4.5 * (x - 4.5) + log(4.5)) * 4;
fplot(y2, [1.5, 9.5])
y3 = @(x) (1/4.5 * (x - 4.5) + log(4.5)) * 4 - 2.5;
fplot(y3, [1.5, 9.5])

scatter([4.5, 2, 2, 2, 8, 8, 8], [y1(4.5), y1(2), y2(2), y3(2), y1(8), y2(8), y3(8)], '.k')
plot([2, 2], [0, y2(2)], '--')
plot([1.75, 2], [y1(2), y1(2)], '-k')
plot([1.75, 2], [y2(2), y2(2)], '-k')
plot([1.75, 2], [y3(2), y3(2)], '-k')
plot([8, 8], [0, y2(8)], '--')
plot([7.75, 8], [y1(8), y1(8)], '-k')
plot([7.75, 8], [y2(8), y2(8)], '-k')
plot([7.75, 8], [y3(8), y3(8)], '-k')

annotation('doublearrow', [1.75 1.75] / 6.65, [y1(2) + 0.4, y2(2) + 0.1] / 10, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [1.75 1.75] / 6.65, [y3(2) + 0.76, y1(2) + 0.4] / 10, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [8 8] / 11.65, [y1(8) - 1.05, y2(8) - 1.25] / 10, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [8 8] / 11.65, [y3(8) - 0.6, y1(8) - 1.05] / 10, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')

text(9.6123, 10.0362, '$l$', 'FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(9.6123, 10.0362 - 2.5, '$l_1$', 'FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(4.6694, 5.6414, '$V_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

text(7.7689, 9.4708, '$V_{B1}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(8.1, 8.0058, '$V_{A1}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(2.1, 2.4387, '$V_{A2}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.7523, 4.1507, '$V_{B2}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

text(1.6, (y1(2) + y2(2))/2, '$\delta$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(7.6, (y1(8) + y2(8))/2, '$\delta$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.6, (y3(2) + y1(2))/2, '$t_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')
text(7.6, (y3(8) + y1(8))/2, '$t_{0}$', 'FontName','Times New Roman', 'Interpreter', 'latex')

xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$V$','FontName','Times New Roman', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'YTick', []);