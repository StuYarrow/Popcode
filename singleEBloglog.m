function h = singleEB(x, y, err, width, label)

le = errorbar(x, y, err);
errorbarlogx(le, width)
set(le, 'linestyle', '-', 'color', 'k', 'marker', 'none')

line([x/(10^(width*0.5)), x, x*(10^(width*0.5))], [y y y], 'color', 'k')

t = text(x * 10^(width * 1.5), y, label);
set(t, 'fontsize', 8, 'verticalalignment', 'middle', 'horizontalalignment', 'left')

h = [le t];