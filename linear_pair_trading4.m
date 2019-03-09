x = [0.5, 2.5];
y = [1, 2 + 1/3];
plot(x, y, 'b')
axis([0,3,0,2.5])
hold on
x = [0.3, 2.5];
y = [0.2, 5 / 3];
plot(x, y, 'b')
scatter([1.5, 1.5, 1.5, 2, 2, 1, 1], [5/3, 1.3, 1, 2, 1 + 1/3, 4/3, 2/3], 'k.')
x = [1.5, 1.5];
y = [5/3, 1];
plot(x, y, '--')

plot([1.45 1.55], [5/3 5/3], '-k')
plot([1.5 1.55], [1.3 1.3], '-k')
plot([1.45 1.5], [1 1], '-k')

plot([1.95 2], [5/3 5/3] + 1/3, '-k')
plot([1.95 2], [1 1] + 1/3, '-k')

plot([0.95 1], [5/3 5/3] - 1/3, '-k')
plot([0.95 1], [1 1] - 1/3, '-k')

plot([1 1], [0 4/3], '--')
plot([2 2], [0 2], '--')
annotation('doublearrow', [0.525 0.525], [1.3 + 5/144, 5/3 - 5/144] / 2.5, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [0.51 0.51], [5/3 - 5/144, 1 + 5/55] / 2.5, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')

annotation('doublearrow', [0.64 0.64], [50/31.3 + 5/144, 1 + 5/55] / 2.5 + 0.11, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [0.38 0.38], [50/31.3 + 5/144, 1 + 5/55] / 2.5 - 0.11, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')

%annotation('doublearrow', [1.475 1.475], [1.3 1])
% xlabel('$P_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
% ylabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')%, 'rotation', 0)
xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$V$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.5, 0.85, '$V_{B0}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.55, 1.2, '$V_{A0}$','FontName','Times New Roman', 'Interpreter', 'latex')

text(2, 1.2 + 0.7, '$V_{A2}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(2, 0.85 + 1/3, '$V_{B2}$','FontName','Times New Roman', 'Interpreter', 'latex')

text(1, 1.2 + 4/3 - 1.3, '$V_{A1}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1, 0.85 - 1/3, '$V_{B1}$','FontName','Times New Roman', 'Interpreter', 'latex')

text(0.4, 1, '$l_{1}$','FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(0.3, 0.35, '$l_{2}$','FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')

text(1.4, (1 + 5/3)/2, '$t_{0}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.55, (1.3 + 5/3)/2, '$\delta$','FontName','Times New Roman', 'Interpreter', 'latex')

text(1.89, (1 + 5/3)/2 + 1/3, '$t_{0}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(0.89, (1 + 5/3)/2 - 1/3, '$t_{0}$','FontName','Times New Roman', 'Interpreter', 'latex')

text(1.89, 0.1, '$P_{B2}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(0.89, 0.1, '$P_{B1}$','FontName','Times New Roman', 'Interpreter', 'latex')
plot([1 1], [0 0.05], '-k', 'LineWidth', 1.3)
plot([2 2], [0 0.05], '-k', 'LineWidth', 1.3)

set(gca, 'XTick', [], 'YTick', []);