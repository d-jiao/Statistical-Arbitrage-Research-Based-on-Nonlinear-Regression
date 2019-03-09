x = [0.5, 2];
y = [1, 2] - (5/3 - 1.3);
plot(x, y, 'b')
axis([0,2.5,0,2.5])
hold on
scatter([1.5 1.5], [5/3 1.3], 'k.')
% xlabel('$P_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
% ylabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')%, 'rotation', 0)
xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$P_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.5, 1.2, '$P$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.4, 1.75, '$P''$','FontName','Times New Roman', 'Interpreter', 'latex')
text(0.5, 0.8, '$l$','FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'YTick', []);