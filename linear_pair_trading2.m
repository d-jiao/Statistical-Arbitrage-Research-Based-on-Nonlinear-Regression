x = [0.5, 2];
y = [1, 2] - (5/3 - 1.3);
plot(x, y, 'b')
axis([0,2.5,0,2.5])
hold on
x = [0.3, 2];
y = [0.2, 4 / 3];
plot(x, y, 'b')
scatter([1.5, 1.5, 1.5], [5/3, 1.3, 1], 'k.')
x = [1.5, 1.5];
y = [5/3, 1];
plot(x, y, '--')
plot([1.45 1.5], [5/3 5/3], '-k')
plot([1.45 1.5], [1.3 1.3], '-k')
plot([1.45 1.5], [1 1], '-k')
annotation('doublearrow', [1.47 1.47] / 2.5, [1.3 + 5/144, 5/3 - 5/144] / 2.5, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
annotation('doublearrow', [1.47 1.47] / 2.5, [1.3 + 5/144, 1 + 5/55] / 2.5, 'Head1Length', 5, 'Head1Width', 5, 'Head2Length', 5, 'Head2Width', 5, 'Units', 'points')
%annotation('doublearrow', [1.475 1.475], [1.3 1])
% xlabel('$P_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
% ylabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')%, 'rotation', 0)
xlabel('$P_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
ylabel('$V$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.5, 1.2, '$V_{A0}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.5, 0.85, '$V_{B}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.4, 1.75, '$V_{A}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(0.5, 0.8, '$l_{1}$','FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(0.3, 0.35, '$l_{2}$','FontSize', 12, 'FontName','Times New Roman', 'Interpreter', 'latex')
text(1.4, (1 + 1.3)/2, '$t_{0}$','FontName','Times New Roman', 'Interpreter', 'latex')
text(1.4, (1.3 + 5/3)/2, '$\delta$','FontName','Times New Roman', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'YTick', []);