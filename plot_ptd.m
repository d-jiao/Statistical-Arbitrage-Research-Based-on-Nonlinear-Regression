for i = 1 : 34
    for j = 1 : i - 1
        scatter(ptd(:,i),ptd(:,j), '.');
        xlabel(universe(i));
        ylabel(universe(j));
        saveas(gcf, [universe{i}, ' & ', universe{j}, '.jpg'], 'jpg');
        [universe{i}, ' & ', universe{j}, '.jpg']
    end
end