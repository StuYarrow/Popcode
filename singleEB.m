function h = singleEB(x, y, err, width, label)

le = errorbar(x, y, err);
errorbarlogx(le, width)
set(le, 'linestyle', '-', 'color', 'k', 'marker', 'none')

line([x/(0.5*10^width), x, x*(0.5 * 10^width)], [y y y], 'color', 'k')

t = text(x * 1.5 * 10^width, y, label);
set(t, 'fontsize', 8, 'verticalalignment', 'middle', 'horizontalalignment', 'left')

h = [le t];