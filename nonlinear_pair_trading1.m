y1 = @(x) exp(x) / 4.5;
fplot(y1, [0.5, 2.5])
hold on
y2 = @(x) (exp(1.5) * (x - 1.5) + exp(1.5))/4.5;
fplot(y2, [0.6, 2.5])
axis([0, 3, 0, 3])
scatter(1.5, exp(1.5) / 4.5, '.k')

text(2.5994, 1.9871, '$l$', 'FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.6657, 0.9638, '$P$', 'FontName','Times New Roman', 'Interpreter', 'latex')

xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$P_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'YTick', []);